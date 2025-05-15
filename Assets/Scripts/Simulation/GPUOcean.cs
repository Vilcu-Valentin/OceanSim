using UnityEngine;
using System.Runtime.InteropServices;
using UnityEngine.Rendering;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class GPUOcean : MonoBehaviour
{
    [Header("Compute Shader & Params")]
    public ComputeShader fftShader;
    public int resolution = 128;
    public float size = 100f;
    public float windSpeed = 10f;
    public float phillipsConstantA = 0.0005f;
    public Vector2 windDirection = new Vector2(1, 0);

    [Header("Scales")]
    public float heightScale = 1f;
    public float chopScale = 1f;

    RenderTexture heightRT, dispRT;
    ComputeBuffer H0_Buffer, Ht_Buffer, Dx_Buffer, Dz_Buffer, TempFFT_H, TempFFT_Dx, TempFFT_Dz, Twiddle_Buffer;

    int initK, updateK, fftHorzK, fftVertK, bitRevK, writeMapsK;
    int log2N, numThreads = 8;

    [StructLayout(LayoutKind.Sequential)]
    struct Complex { public float real, imag; public Complex(float r, float i) { real = r; imag = i; } }

    void Start()
    {
        // -- Mesh Setup --
        var mf = GetComponent<MeshFilter>();
        mf.mesh = GenerateGridMesh();

        // -- Validate & Precompute --
        if (!Mathf.IsPowerOfTwo(resolution)) { enabled = false; return; }
        log2N = (int)Mathf.Log(resolution, 2);

        // -- RenderTextures (bilinear!) --
        heightRT = new RenderTexture(resolution, resolution, 0, RenderTextureFormat.RFloat)
        { enableRandomWrite = true, wrapMode = TextureWrapMode.Repeat, filterMode = FilterMode.Bilinear };
        heightRT.Create();

        dispRT = new RenderTexture(resolution, resolution, 0, RenderTextureFormat.RGFloat)
        { enableRandomWrite = true, wrapMode = TextureWrapMode.Repeat, filterMode = FilterMode.Bilinear };
        dispRT.Create();

        // -- ComputeBuffers --
        int count = resolution * resolution;
        int stride = Marshal.SizeOf(typeof(Complex));
        H0_Buffer = new ComputeBuffer(count, stride);
        Ht_Buffer = new ComputeBuffer(count, stride);
        Dx_Buffer = new ComputeBuffer(count, stride);
        Dz_Buffer = new ComputeBuffer(count, stride);
        TempFFT_H = new ComputeBuffer(count, stride);
        TempFFT_Dx = new ComputeBuffer(count, stride);
        TempFFT_Dz = new ComputeBuffer(count, stride);
        Twiddle_Buffer = new ComputeBuffer(resolution / 2, stride);

        // -- Twiddles --
        CalculateTwiddles(Twiddle_Buffer);

        // -- Kernels --
        initK = fftShader.FindKernel("InitSpectrum");
        updateK = fftShader.FindKernel("UpdateSpectrum");
        fftHorzK = fftShader.FindKernel("FFT_Horizontal_Stage");
        fftVertK = fftShader.FindKernel("FFT_Vertical_Stage");
        fftShader.SetBuffer(fftHorzK, "Twiddles", Twiddle_Buffer);
        fftShader.SetBuffer(fftVertK, "Twiddles", Twiddle_Buffer);
        bitRevK = fftShader.FindKernel("BitReversal_Pass");
        writeMapsK = fftShader.FindKernel("WriteMaps");

        // -- OceanParams CBuffer --
        fftShader.SetInt("_Resolution", resolution);
        fftShader.SetFloat("_Size", size);
        fftShader.SetFloat("_WindSpeed", windSpeed);
        fftShader.SetFloat("_PhillipsA", phillipsConstantA);
        fftShader.SetVector("_WindDir", windDirection.normalized);

        // -- InitSpectrum (H0) --
        fftShader.SetBuffer(initK, "H0", H0_Buffer);
        fftShader.Dispatch(initK, resolution / numThreads, resolution / numThreads, 1);

        // -- Material Setup --
        var mat = GetComponent<MeshRenderer>().material;
        mat.SetTexture("_HeightMap", heightRT);
        mat.SetTexture("_DispMap", dispRT);
        mat.SetFloat("_HeightScale", heightScale);
        mat.SetFloat("_ChopScale", chopScale);
    }

    void Update()
    {
        // -- Update Spectrum Ht --
        fftShader.SetFloat("_Time", Time.time);
        fftShader.SetBuffer(updateK, "H0", H0_Buffer);
        fftShader.SetBuffer(updateK, "Ht", Ht_Buffer);
        fftShader.SetBuffer(updateK, "_DxSpectrum", Dx_Buffer);
        fftShader.SetBuffer(updateK, "_DzSpectrum", Dz_Buffer);
        fftShader.Dispatch(updateK, resolution / numThreads, resolution / numThreads, 1);

        // IFFT the height spectrum into its own buffer:
        var htSpatial = PerformIFFT(Ht_Buffer, TempFFT_H);

        // IFFT the X‐disp spectrum into its own buffer:
        var dxSpatial = PerformIFFT(Dx_Buffer, TempFFT_Dx);

        // IFFT the Z‐disp spectrum into its own buffer:
        var dzSpatial = PerformIFFT(Dz_Buffer, TempFFT_Dz);

        // -- WriteMaps reads SPATIAL Ht_Buffer, Dx_Buffer, Dz_Buffer --
        fftShader.SetBuffer(writeMapsK, "_FFT_Src", htSpatial);       // <— key binding
        fftShader.SetBuffer(writeMapsK, "_DxSpectrum", dxSpatial);
        fftShader.SetBuffer(writeMapsK, "_DzSpectrum", dzSpatial);
        fftShader.SetTexture(writeMapsK, "_HeightRT", heightRT);
        fftShader.SetTexture(writeMapsK, "_DispRT", dispRT);
        fftShader.Dispatch(writeMapsK, resolution / numThreads, resolution / numThreads, 1);
    }

    // Change your signature to return the resulting buffer
    ComputeBuffer PerformIFFT(ComputeBuffer src, ComputeBuffer dst)
    {
        ComputeBuffer a = src, b = dst;

        // -------- horizontal pass (complete 1-D IFFT for every row) --------
        for (int stage = 0; stage < log2N; ++stage)
        {
            fftShader.SetInt("_Stage", stage);
            fftShader.SetInt("_IsInverse", 1);

            fftShader.SetBuffer(fftHorzK, "_FFT_Src", a);
            fftShader.SetBuffer(fftHorzK, "_FFT_Dst", b);
            fftShader.Dispatch(fftHorzK, resolution / numThreads, resolution / numThreads, 1);
            Swap(ref a, ref b);
        }

        // bit-reverse the rows once
        fftShader.SetInt("_BitRevDirection", 0);
        fftShader.SetBuffer(bitRevK, "_FFT_Dst", a);
        fftShader.Dispatch(bitRevK, resolution / numThreads, resolution / numThreads, 1);

        // -------- vertical pass (complete 1-D IFFT for every column) --------
        for (int stage = 0; stage < log2N; ++stage)
        {
            fftShader.SetInt("_Stage", stage);
            fftShader.SetInt("_IsInverse", 1);

            fftShader.SetBuffer(fftVertK, "_FFT_Src", a);
            fftShader.SetBuffer(fftVertK, "_FFT_Dst", b);
            fftShader.Dispatch(fftVertK, resolution / numThreads, resolution / numThreads, 1);
            Swap(ref a, ref b);
        }

        // bit-reverse the columns once
        fftShader.SetInt("_BitRevDirection", 1);
        fftShader.SetBuffer(bitRevK, "_FFT_Dst", a);
        fftShader.Dispatch(bitRevK, resolution / numThreads, resolution / numThreads, 1);

        return a;          // fully spatial-domain, correctly ordered
    }

    void Swap(ref ComputeBuffer x, ref ComputeBuffer y)
    {
        var tmp = x; x = y; y = tmp;
    }

    void CalculateTwiddles(ComputeBuffer buf)
    {
        var data = new Complex[resolution / 2];
        for (int k = 0; k < resolution / 2; k++)
        {
            float ang = -2 * Mathf.PI * k / resolution;
            data[k] = new Complex(Mathf.Cos(ang), Mathf.Sin(ang));
        }
        buf.SetData(data);
    }

    Mesh GenerateGridMesh()
    {
        var mesh = new Mesh { indexFormat = IndexFormat.UInt32 };
        int res = resolution;
        Vector3[] verts = new Vector3[res * res];
        Vector2[] uvs = new Vector2[res * res];
        int[] tris = new int[(res - 1) * (res - 1) * 6];
        float half = size * 0.5f;
        for (int y = 0; y < res; y++) for (int x = 0; x < res; x++)
            {
                int i = y * res + x;
                float u = x / (float)(res - 1), v = y / (float)(res - 1);
                verts[i] = new Vector3(u * size - half, 0, v * size - half);
                uvs[i] = new Vector2(u, v);
            }
        int t = 0;
        for (int y = 0; y < res - 1; y++) for (int x = 0; x < res - 1; x++)
            {
                int i = y * res + x;
                tris[t++] = i; tris[t++] = i + res; tris[t++] = i + res + 1;
                tris[t++] = i; tris[t++] = i + res + 1; tris[t++] = i + 1;
            }
        mesh.vertices = verts; mesh.uv = uvs; mesh.triangles = tris;
        mesh.RecalculateNormals();
        return mesh;
    }

    void OnDestroy()
    {
        heightRT?.Release();
        dispRT?.Release();
        H0_Buffer?.Release();
        Ht_Buffer?.Release();
        Dx_Buffer?.Release();
        Dz_Buffer?.Release();
        Twiddle_Buffer?.Release();
    }
}