using UnityEngine;
using UnityEngine.Formats.Alembic.Importer;

public class AlembicAutoPlayer : MonoBehaviour
{
    // Playback speed multiplier for the animation.
    public float playbackSpeed = 1.0f;

    // Total duration of the animation (in seconds). Set this to your animation's length.
    public float animationDuration = 10.0f;

    // Reference to the AlembicStreamPlayer component.
    private AlembicStreamPlayer streamPlayer;

    // Local time variable to track the current playback time.
    private float currentTime = 0f;

    void Awake()
    {
        // Find the AlembicStreamPlayer component on this GameObject.
        streamPlayer = GetComponent<AlembicStreamPlayer>();
        if (streamPlayer == null)
        {
            Debug.LogError("AlembicStreamPlayer component not found on this GameObject.");
        }
    }

    void Start()
    {
        // Reset the local playback time to 0 and update the stream immediately.
        currentTime = 0f;
        if (streamPlayer != null)
        {
            streamPlayer.UpdateImmediately(currentTime);
        }
    }

    void Update()
    {
        if (streamPlayer != null)
        {
            // Increment the local time based on deltaTime and playbackSpeed.
            currentTime += Time.deltaTime * playbackSpeed;

            // Loop the animation when the currentTime exceeds the animationDuration.
            if (currentTime > animationDuration)
            {
                currentTime -= animationDuration;
            }

            // Update the Alembic animation immediately to the new time.
            streamPlayer.UpdateImmediately(currentTime);
        }
    }
}
