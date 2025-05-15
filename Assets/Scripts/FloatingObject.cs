using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(Rigidbody))]
public class FloatingObject : MonoBehaviour
{
    [Header("Buoyancy Settings")]
    public float offset = 0f;                 // Înălţimea faţă de nivelul apei
    public float floatSpeed = 5f;               // Cât de rapid îşi reglează poziţia
    public float rotateSpeed = 2f;              // Cât de rapid se aliniază la valuri

    [Header("References")]
    public Transform[] floatPoints;             // Cele 4 puncte de testare
    public LayerMask waterLayer;                // Layer-ul apei

    private Vector3[] hitPoints;                // Punctele de impact ale razelor
    private Rigidbody rb;

    void Start()
    {
        rb = GetComponent<Rigidbody>();
        if (floatPoints == null || floatPoints.Length < 3)
            Debug.LogError("Trebuie să ai cel puțin 3 floatPoints configurate!");
        hitPoints = new Vector3[floatPoints.Length];
    }

    void FixedUpdate()
    {
        // 1) Raycast la fiecare punct
        Vector3 avgPoint = Vector3.zero;
        int hits = 0;

        for (int i = 0; i < floatPoints.Length; i++)
        {
            Transform fp = floatPoints[i];
            if (Physics.Raycast(fp.position, Vector3.down, out RaycastHit hit, 100f, waterLayer))
            {
                hitPoints[i] = hit.point;
                avgPoint += hit.point;
                hits++;

                Debug.DrawRay(fp.position, Vector3.down * hit.distance, Color.blue);
            }
            else
            {
                // Dacă nu lovește, păstrează poziția anterioară ca fallback
                hitPoints[i] = fp.position + Vector3.down * (transform.position.y - fp.position.y);
            }
        }

        if (hits == 0) return;

        // 2) Ajustează înălțimea după media punctelor lovite
        avgPoint /= hits;
        Vector3 targetPos = new Vector3(transform.position.x, avgPoint.y + offset, transform.position.z);
        Vector3 moveDir = (targetPos - transform.position) * floatSpeed;
        rb.velocity = new Vector3(rb.velocity.x, moveDir.y, rb.velocity.z);

        // 3) Dacă avem minim 3 hit-uri, calculăm normală de apă și ne aliniem la ea
        if (hits >= 3)
        {
            // Folosim primele 3 puncte pentru plan
            Vector3 v1 = hitPoints[1] - hitPoints[0];
            Vector3 v2 = hitPoints[2] - hitPoints[0];
            Vector3 waterNormal = Vector3.Cross(v1, v2).normalized;

            // Proiecția direcției „forward” pe planul definit de normală
            Vector3 forwardDir = Vector3.ProjectOnPlane(transform.forward, waterNormal).normalized;
            Quaternion targetRot = Quaternion.LookRotation(forwardDir, waterNormal);

            rb.MoveRotation(Quaternion.Slerp(rb.rotation, targetRot, rotateSpeed * Time.fixedDeltaTime));
        }
    }
}
