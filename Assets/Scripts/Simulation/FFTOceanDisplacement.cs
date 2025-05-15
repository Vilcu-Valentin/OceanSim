using UnityEngine;
using System;
using Complex = System.Numerics.Complex;  // alias only Complex to avoid Vector2/3 conflicts
using System.Collections.Generic;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class FFTOceanDisplacement : MonoBehaviour
{
    [Header("Grid & Spectrum")]
    [Tooltip("Must be a power of two.")]
    public int resolution = 64;
    [Tooltip("World-space size of the ocean plane.")]
    public float size = 100f;
    public float windSpeed = 10f;
    [Tooltip("Wind speed (m/s) – only parameter in Phillips spectrum.")]
    public float phillipsConstant = 0.0001f;

    [Header("Amplitude")]
    [Tooltip("Multiplier on all displacements (use 1–10 for testing, then tweak)")]
    public float heightScale = 2f;
    [Tooltip("Multiplier on all displacements (use 1–10 for testing, then tweak)")]
    public float chopScale = 2f;

    [Header("Visualization")]
    public Color deepColor = new Color(0f, 0.1f, 0.3f, 1f);
    public Color peakColor = new Color(1f, 1f, 1f, 1f);

    // --- internals ---
    Mesh mesh;
    Vector3[] baseVerts;
    Vector3[] verts;
    Color[] cols;
    int[] tris;

    // Spectrum storage
    Complex[,] H0;
    float[,] omega;
    Vector2[,] waveK;

    // Working buffers (pre-allocated)
    Complex[,] Ht;
    Complex[,] Hx;
    Complex[,] Hz;

    // FFT twiddle factors: for each stage-length, store twiddles
    Dictionary<int, Complex[]> twiddleCache;

    int N;
    float L;
    const float g = 9.81f;

    void Start()
    {
        // Enforce power of two
        if (!Mathf.IsPowerOfTwo(resolution))
        {
            Debug.LogWarning("Resolution must be a power of two. Rounding to nearest.");
            resolution = Mathf.ClosestPowerOfTwo(resolution);
        }

        N = resolution;
        L = size;

        // Pre-allocate
        H0 = new Complex[N, N];
        omega = new float[N, N];
        waveK = new Vector2[N, N];
        Ht = new Complex[N, N];
        Hx = new Complex[N, N];
        Hz = new Complex[N, N];

        // Build and cache twiddles
        BuildTwiddleCache();

        // Build mesh and spectrum
        BuildMesh();
        InitSpectrum();
    }

    void Update()
    {
        UpdateOcean(Time.time);
    }

    void BuildMesh()
    {
        mesh = new Mesh { name = "FFT Ocean" };

        int vertsPerSide = N + 1;
        baseVerts = new Vector3[vertsPerSide * vertsPerSide];
        verts = new Vector3[baseVerts.Length];
        var uv = new Vector2[baseVerts.Length];
        cols = new Color[baseVerts.Length];
        tris = new int[N * N * 6];

        float half = L * 0.5f;
        float step = L / N;
        int vi = 0;
        for (int z = 0; z <= N; z++)
        {
            for (int x = 0; x <= N; x++, vi++)
            {
                var pos = new Vector3(-half + x * step, 0, -half + z * step);
                baseVerts[vi] = pos;
                verts[vi] = pos;
                uv[vi] = new Vector2((float)x / N, (float)z / N);
                cols[vi] = deepColor;
            }
        }

        int ti = 0;
        for (int z = 0; z < N; z++)
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

        mesh.vertices = verts;
        mesh.uv = uv;
        mesh.triangles = tris;
        mesh.colors = cols;
        mesh.RecalculateBounds(); // keep bounds updated

        GetComponent<MeshFilter>().mesh = mesh;
    }

    void BuildTwiddleCache()
    {
        twiddleCache = new Dictionary<int, Complex[]>();
        int maxLen = N; // FFT length
        for (int len = 2; len <= maxLen; len <<= 1)
        {
            int half = len >> 1;
            var arr = new Complex[half];
            for (int j = 0; j < half; j++)
            {
                double ang = 2 * Math.PI * j / len;
                arr[j] = new Complex(Math.Cos(ang), Math.Sin(ang));
            }
            twiddleCache[len] = arr;
        }
    }

    void InitSpectrum()
    {
        var rand = new System.Random();
        for (int i = 0; i < N; i++)
        {
            float kx = (i <= N / 2 ? i : i - N) * (2f * Mathf.PI / L);
            for (int j = 0; j < N; j++)
            {
                float kz = (j <= N / 2 ? j : j - N) * (2f * Mathf.PI / L);
                var kVec = new Vector2(kx, kz);
                waveK[i, j] = kVec;
                float kLen = kVec.magnitude;

                if (kLen < 1e-6f)
                {
                    H0[i, j] = Complex.Zero;
                    omega[i, j] = 0f;
                }
                else
                {
                    float P = Phillips(kVec, windSpeed, phillipsConstant);
                    double r1 = Gaussian(rand);
                    double r2 = Gaussian(rand);
                    H0[i, j] = new Complex(r1, r2) * Math.Sqrt(P * 0.5);
                    omega[i, j] = Mathf.Sqrt(g * kLen);
                }
            }
        }
        // Enforce Hermitian symmetry
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                int i2 = (N - i) % N;
                int j2 = (N - j) % N;
                H0[i2, j2] = Complex.Conjugate(H0[i, j]);
            }
    }

    void UpdateOcean(float t)
    {
        // 3a) build spectra
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                var h0 = H0[i, j];
                double wt = omega[i, j] * t;
                var expP = Complex.FromPolarCoordinates(1, wt);
                var expM = Complex.FromPolarCoordinates(1, -wt);
                var h = h0 * expP + Complex.Conjugate(h0) * expM;
                Ht[i, j] = h;

                // horizontal displacement
                var k = waveK[i, j]; float kLen = k.magnitude;
                if (kLen > 1e-6f)
                {
                    var factor = Complex.ImaginaryOne * (h / kLen);
                    Hx[i, j] = factor * k.x;
                    Hz[i, j] = factor * k.y;
                }
                else { Hx[i, j] = Complex.Zero; Hz[i, j] = Complex.Zero; }
            }

        // 3b) inverse FFTs (in-place, cached twiddles)
        FFT2D(Ht, true);
        FFT2D(Hx, true);
        FFT2D(Hz, true);

        // 4) apply
        float minH = float.MaxValue, maxH = float.MinValue;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                float v = (float)Ht[i, j].Real;
                minH = Mathf.Min(minH, v);
                maxH = Mathf.Max(maxH, v);
            }

        int idx = 0; float invRange = 1f / (maxH - minH);
        for (int z = 0; z <= N; z++)
            for (int x = 0; x <= N; x++, idx++)
            {
                int i = x % N, j = z % N;
                float h = (float)Ht[i, j].Real * heightScale;
                float dx = (float)Hx[i, j].Real * heightScale * chopScale;
                float dz = (float)Hz[i, j].Real * heightScale * chopScale;

                var bp = baseVerts[idx];
                verts[idx] = new Vector3(bp.x + dx, h, bp.z + dz);

                cols[idx] = Color.Lerp(deepColor, peakColor, (h / heightScale - minH) * invRange);
            }

        mesh.vertices = verts;
        mesh.colors = cols;
        mesh.RecalculateBounds();

        // Compute normals manually
        var norms = new Vector3[verts.Length];
        int side = N + 1;
        for (int z = 0; z <= N; z++)
            for (int x = 0; x <= N; x++)
            {
                int idxN = z * side + x;
                int ixm = Mathf.Max(x - 1, 0), ixp = Mathf.Min(x + 1, N);
                int izm = Mathf.Max(z - 1, 0), izp = Mathf.Min(z + 1, N);
                float hL = verts[z * side + ixm].y;
                float hR = verts[z * side + ixp].y;
                float hD = verts[izm * side + x].y;
                float hU = verts[izp * side + x].y;
                Vector3 n = new Vector3(hL - hR, 2f, hD - hU).normalized;
                norms[idxN] = n;
            }
        mesh.normals = norms;
    }

    // --- FFT helpers using cached twiddles ---
    static int ReverseBits(int n, int bits)
    {
        int rev = 0;
        for (int i = 0; i < bits; i++) { rev = (rev << 1) | (n & 1); n >>= 1; }
        return rev;
    }

    public void FFT(Complex[] buf, bool inverse)
    {
        int n = buf.Length;
        int bits = (int)Math.Log(n, 2);
        for (int i = 0; i < n; i++)
        {
            int j = ReverseBits(i, bits);
            if (j > i) { var tmp = buf[i]; buf[i] = buf[j]; buf[j] = tmp; }
        }

        for (int len = 2; len <= n; len <<= 1)
        {
            int half = len >> 1;
            var tw = twiddleCache[len];
            for (int i = 0; i < n; i += len)
                for (int j = 0; j < half; j++)
                {
                    var u = buf[i + j];
                    var v = buf[i + j + half] * (inverse ? Complex.Conjugate(tw[j]) : tw[j]);
                    buf[i + j] = u + v;
                    buf[i + j + half] = u - v;
                }
        }

        if (inverse)
            for (int i = 0; i < n; i++) buf[i] /= n;
    }

    public void FFT2D(Complex[,] data, bool inverse)
    {
        int n = data.GetLength(0), m = data.GetLength(1);
        var row = new Complex[m];
        for (int i = 0; i < n; i++)
        { for (int j = 0; j < m; j++) row[j] = data[i, j]; FFT(row, inverse); for (int j = 0; j < m; j++) data[i, j] = row[j]; }
        var col = new Complex[n];
        for (int j = 0; j < m; j++)
        { for (int i = 0; i < n; i++) col[i] = data[i, j]; FFT(col, inverse); for (int i = 0; i < n; i++) data[i, j] = col[i]; }
    }

    double Gaussian(System.Random r)
    {
        double u1 = 1 - r.NextDouble(), u2 = 1 - r.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }

    float Phillips(Vector2 k, float windSpd, float A)
    {
        float kLen = k.magnitude;
        if (kLen < 1e-6f) return 0f;
        float k2 = kLen * kLen;
        float L = windSpd * windSpd / g; float L2 = L * L;
        float dotKW = Vector2.Dot(k.normalized, new Vector2(1, 0));
        float damp = 0.0001f; float l2 = (L * damp) * (L * damp);
        return A * Mathf.Exp(-1f / (k2 * L2)) / (k2 * k2) * (dotKW * dotKW) * Mathf.Exp(-k2 * l2);
    }
}
