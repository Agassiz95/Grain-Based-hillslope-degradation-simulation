import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew

# Set the seed for reproducibility
np.random.seed(42)

# Parameters
grid_length = 200
grid_width = 100
origin_x = 50
origin_y = 50
roughness_intensity = 1
num_particles = 1000
horizontal_distance = 100
vertical_distance = 200
stop_percent = 0.7

# Define the slope angle in degrees
slope_angle_degrees = 20  # Adjust this value to control the steepness
# Convert the angle to a slope (rise per unit distance)
slope_y = np.tan(np.radians(slope_angle_degrees))

num_stationary = int(num_particles * 0.15) # How many particles do not move. This is an estimate based on
                                          # field data. In the final histograms there should be ~20% in the zero
                                          # bin for the total histogram and 80% in the inside histogram.
                                          # The average transport distance should be ~6 cm.

# Create the base elevation decreasing along the y-axis
base_elevation = np.arange(grid_length) * -slope_y  # Negative to ensure a downward slope
# Expand base elevation across the x-axis
elevation = np.tile(base_elevation.reshape(-1, 1), (1, grid_width))
# Add surface roughness
roughness = np.random.normal(0, roughness_intensity, (grid_length, grid_width))
elevation += roughness

# Initialize particle positions with the origin at (50, 50)
x_positions = np.full(num_particles, origin_x)
y_positions = np.full(num_particles, origin_y)
paths = [[(origin_x, origin_y)] for _ in range(num_particles)]
has_exited = np.zeros(num_particles, dtype=bool)

# Separate particles into stationary and non-stationary
num_moving = num_particles - num_stationary

stationary_indices = np.arange(num_stationary)
moving_indices = np.arange(num_stationary, num_particles)

# Set up checkpoints relative to the new origin
distance_checkpoints = [i for i in range(-50, 151, 10)]
crossed_particles = {dist: set() for dist in distance_checkpoints}

# Movement directions
moves_downhill = [(0, 1), (-1, 1), (1, 1), (-1, 0), (1, 0)]
moves_radial = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]

# Simulate movement for the moving particles
while np.sum(has_exited) < stop_percent * num_particles:
    move_modes = np.random.choice(['downhill', 'radial'], size=num_moving, p=[0.03, 0.97])

    # Replace Pareto with Exponential distribution
    move_probabilities = np.random.exponential(scale=0.01, size=num_moving)

    move_decisions = move_probabilities < 0.001

    for i, idx in enumerate(moving_indices):
        if has_exited[idx]:
            continue

        if 0 <= y_positions[idx] < grid_length - 1 and move_decisions[i]:
            if move_modes[i] == 'downhill':
                elevation_diffs = []
                for move in moves_downhill:
                    new_x = x_positions[idx] + move[0]
                    new_y = y_positions[idx] + move[1]
                    if 0 <= new_x < grid_width and 0 <= new_y < grid_length:
                        elevation_diffs.append((elevation[y_positions[idx], x_positions[idx]] - elevation[new_y, new_x], move))
                    else:
                        elevation_diffs.append((-np.inf, move))

                elevation_diffs.sort(key=lambda x: x[0], reverse=True)
                chosen_move = elevation_diffs[0][1]
            else:
                # Radial movement bias: favor downward moves
                downhill_bias = 0.05  # Increase to favor downhill movement more
                radial_weights = {
                    (-1, -1): 1 - downhill_bias,  # Up-left
                    (0, -1): 1 - downhill_bias,   # Up
                    (1, -1): 1 - downhill_bias,   # Up-right
                    (-1, 0): 1,  # Left
                    (1, 0): 1,   # Right
                    (-1, 1): 1 + downhill_bias,  # Down-left
                    (0, 1): 1 + downhill_bias,   # Down
                    (1, 1): 1 + downhill_bias    # Down-right
                }

                # Normalize weights
                total_weight = sum(radial_weights.values())
                move_probabilities = [radial_weights[move] / total_weight for move in moves_radial]

                # Select a move based on weighted probabilities
                chosen_move = moves_radial[np.random.choice(len(moves_radial), p=move_probabilities)]

            # Apply the chosen move
            new_x = x_positions[idx] + chosen_move[0]
            new_y = y_positions[idx] + chosen_move[1]

            if 0 <= new_x < grid_width and 0 <= new_y < grid_length:
                x_positions[idx] = new_x
                y_positions[idx] = new_y
                paths[idx].append((new_x, new_y))

            if not (37.5 <= new_x <= 62.5 and 20 <= new_y <= 150):  # x=45,55
                has_exited[idx] = True


# Calculate particle recovery relative to the new origin
for i, path in enumerate(paths):
    for dist in distance_checkpoints:
        grid_dist = origin_y + int((dist / vertical_distance) * grid_length)
        for step_idx in range(1, len(path)):
            prev_y = path[step_idx - 1][1]
            curr_y = path[step_idx][1]
            if prev_y < grid_dist <= curr_y:
                x_pos = path[step_idx][0]
                if origin_x - 12.5 <= x_pos <= origin_x + 12.5:
                    crossed_particles[dist].add(i)
                break

# Print particle counts at checkpoints
print('Particle Counts at Various Distances with New Grid Configuration:')
particle_counts_new_grid = []
for dist, particles in crossed_particles.items():
    count = len(particles)
    particle_counts_new_grid.append(count)
    print(f'At {dist} cm from origin: {count} particles, Recovered: {(count/num_particles)*100:.2f}%')

# Calculate statistics for particle displacement
displacements = []
inside_displacements = []
outside_displacements = []
for path in paths:
    start_x, start_y = path[0]
    end_x, end_y = path[-1]
    displacement = np.sqrt((end_x - start_x) ** 2 + (end_y - start_y) ** 2)
    displacements.append(displacement)
    if 37.5 <= end_x <= 62.5 and 20 <= end_y <= 150: # x=45,55
        inside_displacements.append(displacement)
    else:
        outside_displacements.append(displacement)

displacements = np.array(displacements)
inside_displacements = np.array(inside_displacements)
outside_displacements = np.array(outside_displacements)

mean_displacement = np.mean(displacements)
median_displacement = np.median(displacements)
std_displacement = np.std(displacements)
skew_displacement = skew(displacements)

mean_inside = np.mean(inside_displacements) if inside_displacements.size > 0 else 0
median_inside = np.median(inside_displacements) if inside_displacements.size > 0 else 0
std_inside = np.std(inside_displacements) if inside_displacements.size > 0 else 0
skew_inside = skew(inside_displacements) if inside_displacements.size > 0 else 0

mean_outside = np.mean(outside_displacements) if outside_displacements.size > 0 else 0
median_outside = np.median(outside_displacements) if outside_displacements.size > 0 else 0
std_outside = np.std(outside_displacements) if outside_displacements.size > 0 else 0
skew_outside = skew(outside_displacements) if outside_displacements.size > 0 else 0

print(f"Overall displacement statistics:")
print(f"  Average: {mean_displacement:.2f}")
print(f"  Median: {median_displacement:.2f}")
print(f"  Std Dev: {std_displacement:.2f}")
print(f"  Skewness: {skew_displacement:.2f}")

print(f"\nDisplacement statistics for particles ending inside the boundary:")
print(f"  Average: {mean_inside:.2f}")
print(f"  Median: {median_inside:.2f}")
print(f"  Std Dev: {std_inside:.2f}")
print(f"  Skewness: {skew_inside:.2f}")

print(f"\nDisplacement statistics for particles ending outside the boundary:")
print(f"  Average: {mean_outside:.2f}")
print(f"  Median: {median_outside:.2f}")
print(f"  Std Dev: {std_outside:.2f}")
print(f"  Skewness: {skew_outside:.2f}")

# Report how many grains are still within the sample boundaries
num_inside = len(inside_displacements)
print(f"\nNumber of grains still within the sample boundaries: {num_inside} out of {num_particles}")

plt.figure(figsize=(10, 7))
plt.imshow(elevation, cmap='terrain', origin='lower', extent=[0, grid_width, 0, grid_length])
for path in paths:
    path_x, path_y = zip(*path)
    end_x, end_y = path[-1]
    color = 'red' if 40 <= end_x <= 60 and 30 <= end_y <= 150 else 'black'
    plt.plot(path_x, path_y, color=color, alpha=0.6)

# Add boundary lines
plt.axvline(x=40, color='blue', linestyle='--', label='X Bound Left')
plt.axvline(x=60, color='blue', linestyle='--', label='X Bound Right')
plt.axhline(y=20, color='green', linestyle='-.', label='Y Bound Lower')
plt.axhline(y=150, color='red', linestyle='-.', label='Y Bound Upper')

plt.title('Particle Paths, Origin at (50,50)')
plt.xlabel('X Position (grid units)')
plt.ylabel('Y Position (grid units)')
plt.legend()
plt.colorbar(label='Elevation (cm)')
plt.show()

# Calculate total particles for normalization
total_particles = len(inside_displacements) + len(outside_displacements)

# Calculate total particles for normalization
inside_total = len(inside_displacements)
outside_total = len(outside_displacements)

# Convert frequencies to percentages
bins = [-5, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]
labels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]

# Plot separate histograms for inside and outside, using percentages
plt.figure(figsize=(10, 5))
inside_hist, _ = np.histogram(inside_displacements, bins=bins)
outside_hist, _ = np.histogram(outside_displacements, bins=bins)

# Convert frequencies to percentages based on their own totals
inside_percent = (inside_hist / inside_total) * 100 if inside_total > 0 else np.zeros_like(inside_hist)
outside_percent = (outside_hist / outside_total) * 100 if outside_total > 0 else np.zeros_like(outside_hist)

# Calculate combined histogram
combined_displacements = np.concatenate([inside_displacements, outside_displacements])
combined_hist, _ = np.histogram(combined_displacements, bins=bins)
combined_percent = (combined_hist / total_particles) * 100 if total_particles > 0 else np.zeros_like(combined_hist)

# Plot the combined histogram
plt.figure(figsize=(10, 5))
plt.bar(bins[:-1], combined_percent, width=8, color='blue', alpha=0.6, label='Combined Transport Distances', align='center')
plt.xticks(bins, labels)
plt.xlabel('Transport Distance (cm)')
plt.ylabel('Percentage Recovered (%)')
plt.title('Histogram of Combined Transport Distances (Percentages)')
plt.legend()
plt.show()

# Plot separate histograms for clarity
plt.figure(figsize=(10, 5))
plt.bar(bins[:-1], inside_percent, width=8, color='red', alpha=0.6, label='Inside Bounds', align='edge')
plt.xticks(bins, labels)
plt.xlabel('Transport Distance (cm)')
plt.ylabel('Percentage Recovered (%)')
plt.title('Histogram of Transport Distances - Inside Bounds (Percentages)')
plt.legend()
plt.show()

plt.figure(figsize=(10, 5))
plt.bar(bins[:-1], outside_percent, width=8, color='black', alpha=0.6, label='Outside Bounds', align='edge')
plt.xticks(bins, labels)
plt.xlabel('Transport Distance (cm)')
plt.ylabel('Percentage Recovered (%)')
plt.title('Histogram of Transport Distances - Outside Bounds (Percentages)')
plt.legend()
plt.show()