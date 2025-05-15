using UnityEngine;
using System;
using Complex = System.Numerics.Complex;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;
using System.Collections.Generic;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class FFTOceanDisplacement : MonoBehaviour
{
    [Header("Grid & Spectrum")]
    [Tooltip("Must be a power of two.")]
    public int resolution = 64;
    [Tooltip("World-space size of the ocean plane.")]
    public float size = 100f;
    [Tooltip("Wind speed (m/s) – parameter in Phillips spectrum.")]
    public float windSpeed = 10f;
    [Tooltip("Phillips constant for spectrum tuning.")]
    public float phillipsConstant = 0.0001f;

    [Header("Amplitude")]
    [Tooltip("Height multiplier (1–10).")]
    public float heightScale = 2f;
    [Tooltip("Chop (horizontal) multiplier (1–10).")]
    public float chopScale = 2f;

    [Header("Tiling")]
    [Tooltip("How many tiles in X/Z.")]
    public int tileCount = 3;

    [Header("Visualization")]
    public Color deepColor = new Color(0f, 0.1f, 0.3f, 1f);
    public Color peakColor = new Color(1f, 1f, 1f, 1f);

    float _lastWindSpeed, _lastPhillipsConstant;

    // Internal mesh
    Mesh mesh;
    Vector3[] baseVerts, verts;
    Color[] cols;
    int[] tris;

    // Flattened spectrum arrays
    NativeArray<Complex> H0Arr, HtArr, HxArr, HzArr;
    NativeArray<float> omegaArr;
    NativeArray<Vector2> waveKArr;

    // Cached twiddles for FFT
    Dictionary<int, Complex[]> twiddleCache;

    // Scratch buffer for FFT rows
    Complex[] rowBuf;

    int N;
    float L;
    const float g = 9.81f;

    void Start()
    {
        if (!Mathf.IsPowerOfTwo(resolution))
        {
            Debug.LogWarning("Resolution must be power of two; rounding.");
            resolution = Mathf.ClosestPowerOfTwo(resolution);
        }
        N = resolution;
        L = size;

        int NN = N * N;
        H0Arr = new NativeArray<Complex>(NN, Allocator.Persistent);
        HtArr = new NativeArray<Complex>(NN, Allocator.Persistent);
        HxArr = new NativeArray<Complex>(NN, Allocator.Persistent);
        HzArr = new NativeArray<Complex>(NN, Allocator.Persistent);
        omegaArr = new NativeArray<float>(NN, Allocator.Persistent);
        waveKArr = new NativeArray<Vector2>(NN, Allocator.Persistent);

        rowBuf = new Complex[N];

        BuildTwiddleCache();
        BuildMesh();

        // then spawn the clones
        var root = new GameObject("OceanTiles");
        root.transform.SetParent(transform.parent, false);

        // parent the original under root so it doesn't paint over itself
        transform.SetParent(root.transform, false);

        for (int ix = 0; ix < tileCount; ix++)
            for (int iz = 0; iz < tileCount; iz++)
            {
                if (ix == tileCount / 2 && iz == tileCount / 2) continue; // skip center, that’s you
                var go = Instantiate(gameObject, root.transform);
                go.transform.localPosition = new Vector3(
                  (ix - tileCount / 2) * size,
                  0,
                  (iz - tileCount / 2) * size
                );
                // disable the displacement script on clones so they don’t rebuild the mesh:
                DestroyImmediate(go.GetComponent<FFTOceanDisplacement>());
            }

        _lastWindSpeed = windSpeed;
        _lastPhillipsConstant = phillipsConstant;

        InitSpectrum();
    }

    void OnDestroy()
    {
        if (H0Arr.IsCreated) H0Arr.Dispose();
        if (HtArr.IsCreated) HtArr.Dispose();
        if (HxArr.IsCreated) HxArr.Dispose();
        if (HzArr.IsCreated) HzArr.Dispose();
        if (omegaArr.IsCreated) omegaArr.Dispose();
        if (waveKArr.IsCreated) waveKArr.Dispose();
    }

    void Update()
    {
        // if either value has been tweaked in the Inspector or via script…
        if (!Mathf.Approximately(windSpeed, _lastWindSpeed) ||
            !Mathf.Approximately(phillipsConstant, _lastPhillipsConstant))
        {
            // re-build H0 (and mirror it) using the new values:
            InitSpectrum();

            // update our “last known” so we only do this once per change:
            _lastWindSpeed = windSpeed;
            _lastPhillipsConstant = phillipsConstant;
        }

        // 1) Spectrum update in parallel
        var specJob = new SpectrumJob
        {
            H0 = H0Arr,
            omega = omegaArr,
            waveK = waveKArr,
            Ht = HtArr,
            Hx = HxArr,
            Hz = HzArr,
            time = Time.time,
            chop = chopScale
        };
        JobHandle jh = specJob.Schedule(N * N, 64);
        jh.Complete();

        // 2) Inverse FFT
        FFT2DLinear(HtArr);
        FFT2DLinear(HxArr);
        FFT2DLinear(HzArr);

        // 3) Apply to mesh data
        ApplyToMesh();

        // 4) Upload to mesh
        mesh.SetVertices(verts);
        mesh.SetColors(cols);
        mesh.SetNormals(CalcNormals());
        mesh.RecalculateBounds();
    }

    void BuildMesh()
    {
        mesh = new Mesh { name = "FFT Ocean", indexFormat = UnityEngine.Rendering.IndexFormat.UInt32 };
        int side = N + 1;
        baseVerts = new Vector3[side * side];
        verts = new Vector3[baseVerts.Length];
        var uv = new Vector2[baseVerts.Length];
        cols = new Color[baseVerts.Length];
        tris = new int[N * N * 6];

        float half = L * 0.5f;
        float step = L / N;
        int vi = 0;
        for (int z = 0; z <= N; z++)
            for (int x = 0; x <= N; x++, vi++)
            {
                var p = new Vector3(-half + x * step, 0, -half + z * step);
                baseVerts[vi] = p;
                verts[vi] = p;
                uv[vi] = new Vector2((float)x / N, (float)z / N);
                cols[vi] = deepColor;
            }
        int ti = 0;
        for (int z = 0; z < N; z++)
            for (int x = 0; x < N; x++)
            {
                int i0 = z * side + x;
                tris[ti++] = i0;
                tris[ti++] = i0 + side;
                tris[ti++] = i0 + 1;
                tris[ti++] = i0 + 1;
                tris[ti++] = i0 + side;
                tris[ti++] = i0 + side + 1;
            }
        mesh.vertices = verts;
        mesh.uv = uv;
        mesh.triangles = tris;
        mesh.colors = cols;
        mesh.RecalculateBounds();
        GetComponent<MeshFilter>().sharedMesh = mesh;
    }

    void InitSpectrum()
    {
        var rand = new System.Random();
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                int idx = i * N + j;
                float kx = ((i <= N / 2) ? i : i - N) * (2f * Mathf.PI / L);
                float kz = ((j <= N / 2) ? j : j - N) * (2f * Mathf.PI / L);
                var k = new Vector2(kx, kz);
                waveKArr[idx] = k;
                float kl = k.magnitude;
                if (kl < 1e-6f)
                {
                    H0Arr[idx] = Complex.Zero;
                    omegaArr[idx] = 0f;
                }
                else
                {
                    float P = Phillips(k, windSpeed, phillipsConstant);
                    double r1 = Gaussian(rand), r2 = Gaussian(rand);
                    H0Arr[idx] = new Complex(r1, r2) * Math.Sqrt(P * 0.5);
                    omegaArr[idx] = Mathf.Sqrt(g * kl);
                }
            }
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                int idx = i * N + j;
                int i2 = (N - i) % N, j2 = (N - j) % N;
                H0Arr[i2 * N + j2] = Complex.Conjugate(H0Arr[idx]);
            }
    }

    void FFT2DLinear(NativeArray<Complex> data)
    {
        // Rows
        for (int i = 0; i < N; i++)
        {
            int baseIdx = i * N;
            for (int j = 0; j < N; j++) rowBuf[j] = data[baseIdx + j];
            FFT(rowBuf, true);
            for (int j = 0; j < N; j++) data[baseIdx + j] = rowBuf[j];
        }
        // Columns
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < N; i++) rowBuf[i] = data[i * N + j];
            FFT(rowBuf, true);
            for (int i = 0; i < N; i++) data[i * N + j] = rowBuf[i];
        }
    }

    void ApplyToMesh()
    {
        float minH = float.MaxValue, maxH = float.MinValue;
        for (int idx = 0; idx < N * N; idx++)
        {
            float h = (float)HtArr[idx].Real;
            minH = Mathf.Min(minH, h);
            maxH = Mathf.Max(maxH, h);
        }
        float invR = 1f / (maxH - minH);
        int side = N + 1;
        int v = 0;
        for (int z = 0; z <= N; z++)
            for (int x = 0; x <= N; x++, v++)
            {
                int idx = (x % N) * N + (z % N);
                float h = (float)HtArr[idx].Real * heightScale;
                float dx = (float)HxArr[idx].Real * heightScale * chopScale;
                float dz = (float)HzArr[idx].Real * heightScale * chopScale;
                var bp = baseVerts[v];
                verts[v] = new Vector3(bp.x + dx, h, bp.z + dz);
                float t = ((float)HtArr[idx].Real - minH) * invR;
                cols[v] = Color.Lerp(deepColor, peakColor, t);
            }
    }

    Vector3[] CalcNormals()
    {
        var norms = new Vector3[verts.Length];
        int side = N + 1;
        for (int z = 0; z <= N; z++)
            for (int x = 0; x <= N; x++)
            {
                int i = z * side + x;
                int xm = Mathf.Max(x - 1, 0), xp = Mathf.Min(x + 1, N);
                int zm = Mathf.Max(z - 1, 0), zp = Mathf.Min(z + 1, N);
                float hL = verts[z * side + xm].y;
                float hR = verts[z * side + xp].y;
                float hD = verts[zm * side + x].y;
                float hU = verts[zp * side + x].y;
                norms[i] = new Vector3(hL - hR, 2f, hD - hU).normalized;
            }
        return norms;
    }

    [BurstCompile]
    struct SpectrumJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<Complex> H0;
        [ReadOnly] public NativeArray<float> omega;
        [ReadOnly] public NativeArray<Vector2> waveK;
        public NativeArray<Complex> Ht;
        public NativeArray<Complex> Hx;
        public NativeArray<Complex> Hz;
        public float time;
        public float chop;
        public void Execute(int idx)
        {
            var h0 = H0[idx];
            float w = omega[idx];
            var expP = Complex.FromPolarCoordinates(1, w * time);
            var expM = Complex.FromPolarCoordinates(1, -w * time);
            var h = h0 * expP + Complex.Conjugate(h0) * expM;
            Ht[idx] = h;
            var k = waveK[idx]; float kl = k.magnitude;
            if (kl > 1e-6f)
            {
                var f = Complex.ImaginaryOne * (h / kl);
                Hx[idx] = f * k.x;
                Hz[idx] = f * k.y;
            }
            else { Hx[idx] = Complex.Zero; Hz[idx] = Complex.Zero; }
        }
    }

    static int ReverseBits(int n, int bits)
    {
        int rev = 0;
        for (int i = 0; i < bits; i++) { rev = (rev << 1) | (n & 1); n >>= 1; }
        return rev;
    }

    public void FFT(Complex[] buf, bool inverse)
    {
        int n = buf.Length; int bits = (int)Math.Log(n, 2);
        for (int i = 0; i < n; i++) { int j = ReverseBits(i, bits); if (j > i) { var t = buf[i]; buf[i] = buf[j]; buf[j] = t; } }
        for (int len = 2; len <= n; len <<= 1)
        {
            int half = len >> 1;
            var tw = twiddleCache[len];
            for (int i = 0; i < n; i += len)
                for (int j = 0; j < half; j++)
                {
                    var u = buf[i + j];
                    var v = buf[i + j + half] * (inverse ? Complex.Conjugate(tw[j]) : tw[j]);
                    buf[i + j] = u + v; buf[i + j + half] = u - v;
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
        {
            for (int j = 0; j < m; j++) row[j] = data[i, j]; FFT(row, inverse); for (int j = 0; j < m; j++) data[i, j] = row[j];
        }
        var col = new Complex[n];
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++) col[i] = data[i, j]; FFT(col, inverse); for (int i = 0; i < n; i++) data[i, j] = col[i];
        }
    }

    void BuildTwiddleCache()
    {
        twiddleCache = new Dictionary<int, Complex[]>();
        for (int len = 2; len <= N; len <<= 1)
        {
            int half = len >> 1;
            var arr = new Complex[half];
            for (int j = 0; j < half; j++) { double ang = 2 * Math.PI * j / len; arr[j] = new Complex(Math.Cos(ang), Math.Sin(ang)); }
            twiddleCache[len] = arr;
        }
    }

    double Gaussian(System.Random r)
    {
        double u1 = 1 - r.NextDouble(), u2 = 1 - r.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }

    float Phillips(Vector2 k, float wind, float A)
    {
        float kl = k.magnitude; if (kl < 1e-6f) return 0f;
        float k2 = kl * kl;
        float Lw = wind * wind / g, L2 = Lw * Lw;
        float dk = Vector2.Dot(k.normalized, new Vector2(1, 0));
        float damp = 0.0001f, l2 = (Lw * damp) * (Lw * damp);
        return A * Mathf.Exp(-1f / (k2 * L2)) / (k2 * k2) * (dk * dk) * Mathf.Exp(-k2 * l2);
    }

    /// <summary>
    /// Samples the ocean height (world-space Y) at an arbitrary world XZ position,
    /// by bilinear-interpolating your HtArr (post-FFT) grid.
    /// </summary>
    public float SampleHeight(Vector3 worldPos)
    {
        // 1) Transform into local (mesh) space
        Vector3 local = transform.InverseTransformPoint(worldPos);
        // 2) Map X/Z from [-L/2, +L/2] to [0, N)
        float u = (local.x + L * 0.5f) / L * N;
        float v = (local.z + L * 0.5f) / L * N;
        // 3) Get integer cell and frac
        int i0 = Mathf.FloorToInt(u) % N;
        int j0 = Mathf.FloorToInt(v) % N;
        if (i0 < 0) i0 += N;
        if (j0 < 0) j0 += N;
        int i1 = (i0 + 1) % N;
        int j1 = (j0 + 1) % N;
        float fu = u - Mathf.Floor(u);
        float fv = v - Mathf.Floor(v);

        // 4) Fetch four corner heights
        float h00 = (float)HtArr[i0 * N + j0].Real * heightScale;
        float h10 = (float)HtArr[i1 * N + j0].Real * heightScale;
        float h01 = (float)HtArr[i0 * N + j1].Real * heightScale;
        float h11 = (float)HtArr[i1 * N + j1].Real * heightScale;

        // 5) Bilinear interp
        float h0 = Mathf.Lerp(h00, h10, fu);
        float h1 = Mathf.Lerp(h01, h11, fu);
        return Mathf.Lerp(h0, h1, fv) + transform.position.y;
    }

    /// <summary>
    /// Samples the ocean normal at an arbitrary world XZ position,
    /// by bilinear-interpolating your HxArr/HzArr (slopes) and constructing the normal.
    /// </summary>
    public Vector3 SampleNormal(Vector3 worldPos)
    {
        // same UV calc as above
        Vector3 local = transform.InverseTransformPoint(worldPos);
        float u = (local.x + L * 0.5f) / L * N;
        float v = (local.z + L * 0.5f) / L * N;
        int i0 = Mathf.FloorToInt(u) % N;
        int j0 = Mathf.FloorToInt(v) % N;
        if (i0 < 0) i0 += N;
        if (j0 < 0) j0 += N;
        int i1 = (i0 + 1) % N;
        int j1 = (j0 + 1) % N;
        float fu = u - Mathf.Floor(u);
        float fv = v - Mathf.Floor(v);

        // fetch dx, dz slopes
        float dx00 = (float)HxArr[i0 * N + j0].Real * heightScale * chopScale;
        float dx10 = (float)HxArr[i1 * N + j0].Real * heightScale * chopScale;
        float dx01 = (float)HxArr[i0 * N + j1].Real * heightScale * chopScale;
        float dx11 = (float)HxArr[i1 * N + j1].Real * heightScale * chopScale;

        float dz00 = (float)HzArr[i0 * N + j0].Real * heightScale * chopScale;
        float dz10 = (float)HzArr[i1 * N + j0].Real * heightScale * chopScale;
        float dz01 = (float)HzArr[i0 * N + j1].Real * heightScale * chopScale;
        float dz11 = (float)HzArr[i1 * N + j1].Real * heightScale * chopScale;

        // bilinear interp on slopes
        float dx0 = Mathf.Lerp(dx00, dx10, fu);
        float dx1 = Mathf.Lerp(dx01, dx11, fu);
        float dz0 = Mathf.Lerp(dz00, dz10, fu);
        float dz1 = Mathf.Lerp(dz01, dz11, fu);

        float dx = Mathf.Lerp(dx0, dx1, fv);
        float dz = Mathf.Lerp(dz0, dz1, fv);

        // normal from slope
        Vector3 nLocal = new Vector3(-dx, 1f, -dz).normalized;
        return transform.TransformDirection(nLocal);
    }

}
