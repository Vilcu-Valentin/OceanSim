using UnityEngine;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class PhillipsGerstnerOcean : MonoBehaviour
{
    [Header("Grid Settings")]
    [Tooltip("Number of quads per side (must be ≥ 1)")]
    public int resolution = 64;
    [Tooltip("Size of the ocean plane in world units")]
    public float size = 100f;

    [Header("Spectral Settings")]
    [Tooltip("Wind speed (m/s) controls overall wave height")]
    public float windSpeed = 12f;
    [Tooltip("Direction the wind blows (unit vector)")]
    public Vector2 windDirection = new Vector2(1, 0);
    [Tooltip("Phillips spectrum constant (small ~1e-3…1e-1)")]
    public float phillipsConstant = 0.0005f;

    [Header("Wave Composition")]
    [Tooltip("Number of Gerstner waves to sum")]
    [Range(1, 32)]
    public int waveCount = 16;
    [Tooltip("Min and max wavelength for sampling")]
    public float minWavelength = 5f, maxWavelength = 50f;

    // internal
    struct Wave { public Vector2 dir; public float k; public float A; public float omega; public float phase; }
    Wave[] waves;
    Mesh mesh;
    Vector3[] baseVerts;

    void Start()
    {
        BuildMesh();
        InitWaves();
    }

    void Update()
    {
        UpdateWaves();
    }

    // 1) Build a flat grid mesh, store its base vertices
    void BuildMesh()
    {
        mesh = new Mesh();
        mesh.name = "OceanMesh";
        int vertsPerSide = resolution + 1;
        baseVerts = new Vector3[vertsPerSide * vertsPerSide];
        Vector2[] uv = new Vector2[baseVerts.Length];
        int[] tris = new int[resolution * resolution * 6];

        float half = size * 0.5f;
        float step = size / resolution;
        int vi = 0;
        for (int z = 0; z <= resolution; z++)
        {
            for (int x = 0; x <= resolution; x++)
            {
                baseVerts[vi] = new Vector3(-half + x * step, 0f, -half + z * step);
                uv[vi] = new Vector2((float)x / resolution, (float)z / resolution);
                vi++;
            }
        }

        int ti = 0;
        for (int z = 0; z < resolution; z++)
        {
            for (int x = 0; x < resolution; x++)
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
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();

        var mf = GetComponent<MeshFilter>();
        mf.mesh = mesh;
    }

    // 2) Initialize wave components by sampling the Phillips spectrum
    void InitWaves()
    {
        waves = new Wave[waveCount];
        Vector2 windDirN = windDirection.normalized;
        float g = 9.81f;
        float L = windSpeed * windSpeed / g;  // largest wave scale

        for (int i = 0; i < waveCount; i++)
        {
            // random direction & wavelength in [min,max]
            float theta = Random.Range(0f, Mathf.PI * 2f);
            Vector2 dir = new Vector2(Mathf.Cos(theta), Mathf.Sin(theta));

            float λ = Random.Range(minWavelength, maxWavelength);
            float k = 2 * Mathf.PI / λ;

            // Phillips spectrum P(k)
            float kd = Vector2.Dot(dir, windDirN);
            kd = Mathf.Max(kd, 0f);
            float phillips = phillipsConstant
                            * Mathf.Exp(-1f / (k * L * k * L))
                            / (k * k * k * k)
                            * (kd * kd);

            float A = Mathf.Sqrt(phillips * 0.5f);
            float ω = Mathf.Sqrt(g * k);
            float φ = Random.Range(0f, Mathf.PI * 2f);

            waves[i] = new Wave { dir = dir, k = k, A = A, omega = ω, phase = φ };
        }
    }

    // 3) Each frame, displace mesh verts (vertical + horizontal) and recalc normals
    void UpdateWaves()
    {
        var verts = new Vector3[baseVerts.Length];
        baseVerts.CopyTo(verts, 0);
        float t = Time.time;

        for (int i = 0; i < verts.Length; i++)
        {
            Vector3 p = baseVerts[i];
            Vector2 xz = new Vector2(p.x, p.z);
            Vector3 disp = Vector3.zero;

            // sum Gerstner contributions
            foreach (var w in waves)
            {
                float d = Vector2.Dot(w.dir, xz);
                float phase = d * w.k - w.omega * t + w.phase;
                float cos = Mathf.Cos(phase);
                float sin = Mathf.Sin(phase);

                // horizontal displacement
                disp.x += -w.dir.x * w.A * cos;
                disp.z += -w.dir.y * w.A * cos;
                // vertical displacement
                disp.y += w.A * sin;
            }

            verts[i] = p + disp;
        }

        mesh.vertices = verts;
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
    }
}
