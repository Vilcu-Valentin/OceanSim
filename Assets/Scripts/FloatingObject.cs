using UnityEngine;

[RequireComponent(typeof(Rigidbody))]
public class FloatingObject : MonoBehaviour
{
    [Header("Buoyancy Settings")]
    public float offset = 0f;    // height above wave
    public float floatSpeed = 5f;
    public float rotateSpeed = 2f;

    [Header("References")]
    public Transform[] floatPoints;          // test points on your hull
    public FFTOceanDisplacement oceanSolver; // drag in your FFT Ocean GameObject

    private Vector3[] samplePositions;
    private Vector3[] sampleNormals;
    private Rigidbody rb;

    void Start()
    {
        rb = GetComponent<Rigidbody>();
        if (floatPoints == null || floatPoints.Length < 3)
            Debug.LogError("Need at least 3 floatPoints!");
        samplePositions = new Vector3[floatPoints.Length];
        sampleNormals = new Vector3[floatPoints.Length];
    }

    void FixedUpdate()
    {
        int hits = floatPoints.Length;
        Vector3 avgPos = Vector3.zero;
        Vector3 avgNormal = Vector3.zero;

        // 1) sample height & normal at each floatPoint
        for (int i = 0; i < floatPoints.Length; i++)
        {
            var fp = floatPoints[i].position;
            float h = oceanSolver.SampleHeight(fp);
            Vector3 n = oceanSolver.SampleNormal(fp);
            samplePositions[i] = new Vector3(fp.x, h, fp.z);
            sampleNormals[i] = n;

            avgPos += samplePositions[i];
            avgNormal += n;
        }

        // 2) average
        avgPos /= hits;
        avgNormal = (avgNormal / hits).normalized;

        // 3) move rigidbody to match wave height
        Vector3 targetPos = new Vector3(transform.position.x,
                                        avgPos.y + offset,
                                        transform.position.z);
        Vector3 vel = rb.velocity;
        vel.y = (targetPos.y - transform.position.y) * floatSpeed;
        rb.velocity = vel;

        // 4) align rotation to average normal
        Vector3 forwardOnPlane = Vector3.ProjectOnPlane(transform.forward, avgNormal);
        Quaternion targetRot = Quaternion.LookRotation(forwardOnPlane, avgNormal);
        rb.MoveRotation(Quaternion.Slerp(rb.rotation,
                                         targetRot,
                                         rotateSpeed * Time.fixedDeltaTime));
    }
}
