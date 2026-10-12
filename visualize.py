import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Load image
image = np.loadtxt("/Users/jonahchen/physics427-2025-project-1-the-best-team/bh_image.pgm", skiprows=3)

M = 4
ny, nx = image.shape
alpha = np.linspace(-nx/(2*M), nx/(2*M), nx)
beta = np.linspace(-ny/(2*M), ny/(2*M), ny)

# --- Map each pixel independently ---
# Zeros stay black, non-zero pixels get mapped red->orange->white
stretched_image = np.copy(image)

# Avoid global min/max scaling; scale non-zero pixels relative to the max only
max_val = np.max(stretched_image)
stretched_image[stretched_image > 0] = stretched_image[stretched_image > 0] / max_val

# --- Custom colormap ---
colors = [(0, "black"), (0.001, "red"), (0.5, "orange"), (1, "white")]
custom_cmap = LinearSegmentedColormap.from_list("black_red_orange_white", colors)

plt.imshow(stretched_image, cmap=custom_cmap, origin='lower',
           extent=[alpha[0], alpha[-1], beta[0], beta[-1]])
plt.colorbar(label="Intensity (0 is black)")
plt.title("Black Hole Raytrace, spin = 0.3")
plt.xlabel("α / M")
plt.ylabel("β / M")
plt.show()
