using UnityEngine;
using System;
using System.Numerics;  // for Complex

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class FFTOceanDisplacement : MonoBehaviour
{
    [Header("Grid & Spectrum")]
    [Tooltip("Must be a power of two.")]
    public int resolution = 64;
    [Tooltip("World‐space size of the ocean plane.")]
    public float size = 100f;
    [Tooltip("Wind speed (m/s) – only parameter in Phillips spectrum.")]
    public float windSpeed = 10f;
    [Tooltip("Phillips constant (tune for larger/smaller waves).")]
    public float phillipsConstant = 0.0001f;

    [Header("Amplitude")]
    [Tooltip("Multiplier on all displacements (use 1–10 for testing, then tweak)")]
    public float heightScale = 2f;
    [Tooltip("Multiplier on all displacements (use 1–10 for testing, then tweak)")]
    public float chopScale = 2f;

    [Header("Visualization")]
    [Tooltip("Color at wave troughs.")]
    public Color deepColor = new Color(0f, 0.1f, 0.3f, 1f);
    [Tooltip("Color at wave peaks.")]
    public Color peakColor = new Color(1f, 1f, 1f, 1f);

    // --- internals ---
    Mesh mesh;
    UnityEngine.Vector3[] baseVerts;
    Complex[,] H0;        // initial spectrum
    float[,] omega;      // angular frequencies ω(k)
    UnityEngine.Vector2[,] waveK;     // k vectors
    int N;                // alias for resolution
    float L;              // alias for size
    const float g = 9.81f;

    void Start()
    {
        // enforce power of two
        if (!Mathf.IsPowerOfTwo(resolution))
        {
            Debug.LogWarning("Resolution must be a power of two. Rounding to nearest.");
            resolution = Mathf.ClosestPowerOfTwo(resolution);
        }

        N = resolution;
        L = size;

        BuildMesh();
        InitSpectrum();
    }

    void Update()
    {
        UpdateOcean(Time.time);
    }

    // 1) Build a flat (N+1)x(N+1) grid and assign initial deepColor
    void BuildMesh()
    {
        mesh = new Mesh();
        mesh.name = "FFT Ocean";

        int vertsPerSide = N + 1;
        baseVerts = new UnityEngine.Vector3[vertsPerSide * vertsPerSide];
        UnityEngine.Vector2[] uv = new UnityEngine.Vector2[baseVerts.Length];
        Color[] cols = new Color[baseVerts.Length];
        int[] tris = new int[N * N * 6];

        float half = L * 0.5f;
        float step = L / N;
        int vi = 0;
        for (int z = 0; z <= N; z++)
        {
            for (int x = 0; x <= N; x++, vi++)
            {
                baseVerts[vi] = new UnityEngine.Vector3(-half + x * step, 0f, -half + z * step);
                uv[vi] = new UnityEngine.Vector2((float)x / N, (float)z / N);
                cols[vi] = deepColor;
            }
        }

        int ti = 0;
        for (int z = 0; z < N; z++)
        {
            for (int x = 0; x < N; x++)
            {
                int i0 = z * vertsPerSide + x;
                tris[ti++] = i0;
                tris[ti++] = i0 + vertsPerSide;
                tris[ti++] = i0 + 1;

                tris[ti++] = i0 + 1;
                tris[ti++] = i0 + vertsPerSide;
                tris[ti++] = i0 + vertsPerSide + 1;
            }
        }

        mesh.vertices = baseVerts;
        mesh.uv = uv;
        mesh.triangles = tris;
        mesh.colors = cols;
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();

        var mf = GetComponent<MeshFilter>();
        mf.mesh = mesh;
    }

    // 2) Initialize H0(k) and ω(k) for each grid cell, then enforce Hermitian symmetry
    void InitSpectrum()
    {
        H0 = new Complex[N, N];
        omega = new float[N, N];
        waveK = new UnityEngine.Vector2[N, N];
        var rand = new System.Random();

        // build spectrum
        for (int i = 0; i < N; i++)
        {
            float kx = (i <= N / 2 ? i : i - N) * (2f * Mathf.PI / L);
            for (int j = 0; j < N; j++)
            {
                float kz = (j <= N / 2 ? j : j - N) * (2f * Mathf.PI / L);
                var kVec = new UnityEngine.Vector2(kx, kz);
                waveK[i, j] = kVec;
                float kLen = kVec.magnitude;

                if (kLen < 1e-6f)
                {
                    H0[i, j] = Complex.Zero;
                    omega[i, j] = 0f;
                }
                else
                {
                    // Phillips spectrum
                    float P = Phillips(kVec, windSpeed, phillipsConstant);
                    // Gaussian random pair
                    double r1 = Gaussian(rand);
                    double r2 = Gaussian(rand);
                    // H0 = (ξ_r + i ξ_i) * sqrt(P/2)
                    H0[i, j] = new Complex(r1, r2) * Math.Sqrt(P * 0.5);
                    omega[i, j] = Mathf.Sqrt(g * kLen);
                }
            }
        }

        // Hermitian symmetry: H0[-k] = conj(H0[k]) to ensure real h(x)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                int i2 = (N - i) % N;
                int j2 = (N - j) % N;
                H0[i2, j2] = Complex.Conjugate(H0[i, j]);
            }
    }

    // Box–Muller for a single N(0,1) sample
    double Gaussian(System.Random r)
    {
        double u1 = 1.0 - r.NextDouble();
        double u2 = 1.0 - r.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // Phillips spectrum with a little damping to tame very short waves
    float Phillips(UnityEngine.Vector2 k, float windSpd, float A)
    {
        float kLen = k.magnitude;
        float k2 = kLen * kLen;
        float L = windSpd * windSpd / g;
        float L2 = L * L;
        float dotKW = UnityEngine.Vector2.Dot(k.normalized, new UnityEngine.Vector2(1, 0)); // assume wind along +X
        float damp = 0.0001f;
        float l2 = L * damp * L * damp;

        if (kLen < 1e-6f) return 0f;

        float ph = A * Mathf.Exp(-1f / (k2 * L2)) / (k2 * k2)
                   * (dotKW * dotKW)
                   * Mathf.Exp(-k2 * l2);
        return ph;
    }

    // 3) Each frame: build H(t), then do 2D inverse FFTs for height & horizontal displacements
    void UpdateOcean(float t)
    {
        // temporary spectra
        var Ht = new Complex[N, N];
        var Hx = new Complex[N, N];
        var Hz = new Complex[N, N];

        // 3a) H˜(k, t) = H0 e^{iωt} + conj(H0) e^{-iωt}
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                var h0 = H0[i, j];
                float w = omega[i, j];
                var expP = Complex.Exp(new Complex(0, w * t));
                var expM = Complex.Exp(new Complex(0, -w * t));
                var h = h0 * expP + Complex.Conjugate(h0) * expM;
                Ht[i, j] = h;

                // horizontal: i*(k/|k|)*h
                var k = waveK[i, j];
                float kLen = k.magnitude;
                if (kLen > 1e-6f)
                {
                    var factor = new Complex(0, 1) * (h / kLen);
                    Hx[i, j] = factor * k.x;
                    Hz[i, j] = factor * k.y;
                }
            }

        // 3b) inverse FFT → real height & displacement fields
        FFT2D(Ht, true);
        FFT2D(Hx, true);
        FFT2D(Hz, true);

        // 4) apply to mesh, recalc normals & color by height
        var verts = mesh.vertices;
        var cols = mesh.colors;

        // find min/max height for vertex coloring
        float minH = float.MaxValue, maxH = float.MinValue;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                float v = (float)Ht[i, j].Real;
                if (v < minH) minH = v;
                if (v > maxH) maxH = v;
            }

        Debug.Log($"[FFT Ocean] raw height range = {minH:F4} → {maxH:F4}; scaled = {minH * heightScale:F2} → {maxH * heightScale:F2}");

        int idx = 0;
        for (int z = 0; z <= N; z++)
        {
            for (int x = 0; x <= N; x++, idx++)
            {
                int i = x % N;
                int j = z % N;

                // amplify everything so it's visible
                float h = (float)Ht[i, j].Real * heightScale;
                float dx = (float)Hx[i, j].Real * heightScale * chopScale;
                float dz = (float)Hz[i, j].Real * heightScale * chopScale;

                var bp = baseVerts[idx];
                verts[idx] = new UnityEngine.Vector3(bp.x + dx, h, bp.z + dz);

                float tcol = Mathf.InverseLerp(minH, maxH, h);
                cols[idx] = Color.Lerp(deepColor, peakColor, tcol);
            }
        }

        mesh.vertices = verts;
        mesh.colors = cols;
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
    }

    // --- FFT helpers ---

    static int ReverseBits(int n, int bits)
    {
        int rev = 0;
        for (int i = 0; i < bits; i++)
        {
            rev = (rev << 1) | (n & 1);
            n >>= 1;
        }
        return rev;
    }

    // 1D Cooley–Tuk FFT (in-place). inverse = true → IFFT
    public static void FFT(Complex[] buffer, bool inverse)
    {
        int n = buffer.Length;
        int bits = (int)Math.Log(n, 2);

        // bit reversal
        for (int i = 0; i < n; i++)
        {
            int j = ReverseBits(i, bits);
            if (j > i)
            {
                var tmp = buffer[i];
                buffer[i] = buffer[j];
                buffer[j] = tmp;
            }
        }

        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = 2 * Math.PI / len * (inverse ? -1 : 1);
            var wlen = new Complex(Math.Cos(ang), Math.Sin(ang));

            for (int i = 0; i < n; i += len)
            {
                var w = Complex.One;
                int half = len >> 1;
                for (int j = 0; j < half; j++)
                {
                    var u = buffer[i + j];
                    var v = buffer[i + j + half] * w;
                    buffer[i + j] = u + v;
                    buffer[i + j + half] = u - v;
                    w *= wlen;
                }
            }
        }

        if (inverse)
        {
            for (int i = 0; i < n; i++)
                buffer[i] /= n;
        }
    }

    // 2D FFT via row-wise then column-wise 1D FFT
    public static void FFT2D(Complex[,] data, bool inverse)
    {
        int n = data.GetLength(0);
        int m = data.GetLength(1);

        var row = new Complex[m];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++) row[j] = data[i, j];
            FFT(row, inverse);
            for (int j = 0; j < m; j++) data[i, j] = row[j];
        }

        var col = new Complex[n];
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++) col[i] = data[i, j];
            FFT(col, inverse);
            for (int i = 0; i < n; i++) data[i, j] = col[i];
        }
    }
}
