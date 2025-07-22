import math
import klayout.db as pya

# Parameters in micrometers (um)
L = 22.5         # arc length in um
chord = 5.0      # distance between the endpoints' centers in um
w = 0.45         # path width in um

# Solve for theta using bisection:
# Given: chord = 2*(L/theta)*sin(theta/2)
def solve_theta(L, chord):
    low = 0.001
    high = math.pi  # reasonable upper bound since chord <= L
    for _ in range(100):
        theta = (low + high) / 2
        f = 2 * (L/theta) * math.sin(theta/2) - chord
        if f > 0:
            low = theta
        else:
            high = theta
    return (low + high) / 2

theta = solve_theta(L, chord)
r = L/theta  # radius in um

print(f"theta = {theta} rad, r = {r} um")

# Create an arc (DPath) in KLayout. The arc will have the calculated parameters.
# We generate many points along the arc.
num_points = 100
points = []

# Determine the start angle so that the chord is horizontal.
# Here we choose the arc so that its chord lies along the x-axis.
# The center of the circle is at (0, r) in um.
start_angle = -math.pi/2 - theta/2
for i in range(num_points + 1):
    angle = start_angle + i * (theta/num_points)
    x = r * math.cos(angle)
    y = r * math.sin(angle) + r  # shift circle center to (0, r)
    points.append(pya.DPoint(x, y))  # points in um

# 計算雖然由參數生成的弧線實際長度
arc_length = 0
for i in range(1, len(points)):
    dx = points[i].x - points[i-1].x
    dy = points[i].y - points[i-1].y
    arc_length += math.sqrt(dx*dx + dy*dy)
print("原始弧線長度 (um):", arc_length)

# 假設我們希望曲線長度為一個整數（取最接近的整數）
target_length = round(arc_length)
print("目標弧線長度 (um):", target_length)

# 計算調整縮放因子
scale = target_length / arc_length
print("縮放因子:", scale)

# 用此縮放因子調整所有點的座標
points = [ pya.DPoint(pt.x * scale, pt.y * scale) for pt in points]

# 驗證新弧線長度是否為整數
new_length = 0
for i in range(1, len(points)):
    dx = points[i].x - points[i-1].x
    dy = points[i].y - points[i-1].y
    new_length += math.sqrt(dx*dx + dy*dy)
print("調整後弧線長度 (um):", new_length)

# Create a DPath from these points with width 'w' (in um)
arc_path = pya.DPath(points, w)

# Convert to database units (nm):
layout = pya.Layout()
layout.dbu = 0.001  # 1 dbu = 0.001 um ==> 1 um = 1000 nm
top_cell = layout.create_cell("TOP")
layer_index = layout.layer(1, 0)
top_cell.shapes(layer_index).insert(arc_path)

# Save layout to GDS
layout.write("output.gds")