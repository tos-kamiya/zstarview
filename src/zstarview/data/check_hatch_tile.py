"""
Generates a debug image for a diagonal hatch pattern.

This script is a simple utility to create and save a grayscale hatch tile
as a PNG image. It's useful for visually inspecting the pattern generated
by the hatching algorithm used elsewhere in the application.
"""

import math

import numpy as np
from PIL import Image


def make_hatch_tile_pil(W: int, H: int, line_px: int, strength: int) -> Image.Image:
    """
    Generates a hatch tile image and returns it as a PIL.Image (grayscale).

    This function creates a grayscale image with a diagonal line pattern (hatch).
    The lines are anti-aliased by calculating the distance of each pixel
    from a central diagonal line.

    Args:
        W: The width of the tile in pixels.
        H: The height of the tile in pixels.
        line_px: The thickness of the hatch lines in pixels.
        strength: The grayscale value (0-255) of the hatch lines.

    Returns:
        A Pillow `Image` object in grayscale ("L") mode.
    """
    # --- Calculate line thickness parameters ---
    # Normalize line thickness based on tile dimensions to maintain appearance.
    norm = math.sqrt(W * W + H * H)
    P = W * H
    band_u = max(1, int(round(line_px * norm)))

    # --- Create a coordinate grid ---
    xs = np.arange(W, dtype=np.int32)[None, :]
    ys = np.arange(H, dtype=np.int32)[:, None]

    # --- Generate the diagonal pattern ---
    # Calculate the distance of each pixel from the main diagonal.
    u = H * xs - W * ys
    u_mod = np.mod(u, P)
    dist = np.minimum(u_mod, P - u_mod)
    # Create a boolean mask for pixels that fall within the line band.
    mask = dist <= (band_u / 2)

    # --- Create the image array ---
    # Start with a black image.
    arr = np.zeros((H, W), dtype=np.uint8)
    # Set the color of the pixels under the mask to the desired strength.
    arr[mask] = np.uint8(np.clip(strength, 0, 255))

    return Image.fromarray(arr, mode="L")


# --- Main execution: Generate and save a sample tile ---
if __name__ == "__main__":
    # Create a sample 32x20 tile with a line thickness of 5px and a color of 200.
    tile = make_hatch_tile_pil(32, 20, 5, 200)
    # Save the generated tile to a file for visual inspection.
    output_filename = "hatch_tile_debug.png"
    tile.save(output_filename)
    print(f"Saved sample hatch tile to '{output_filename}'")
