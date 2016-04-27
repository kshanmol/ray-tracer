from collections import namedtuple
from math import log10

Point = namedtuple('Point', ['x', 'y', 'z'])
Plane = namedtuple('Plane', ['A', 'B', 'C', 'D'])
Square = namedtuple('Square', ['A', 'B', 'C', 'D'])

SQUARES_PER_UNIT = 10
"""

       A _______________________ B
        /                      /
       /                      /
      /                      /
     /                      /
    /                      /
   /______________________/
  D                       C
"""
round_amount = int(log10(SQUARES_PER_UNIT))
step = 1.0/SQUARES_PER_UNIT

def float_stepper(from_, to_, step_):
    while from_ <= to_:
        from_ = round(from_, round_amount)
        yield from_
        from_ += step_

plane = Plane(A = Point(-5, -1, -5),    # A
              B = Point(5, -1, -5),     # B
              C = Point(5, -1, 5),      # C
              D = Point(-5, -1, 5),     # D
)

slices_x = [plane.A.x]
squares = []
for slice_x in float_stepper(plane.A.x + step, plane.B.x + step, step):
    slices_x.append(slice_x)
    slices_z = [plane.A.z]
    for slice_z in float_stepper(plane.A.z + step, plane.D.z + step, step):
        slices_z.append(slice_z)
        squares.append(Square(
            A=Point(slices_x[-2], -1, slices_z[-2]),
            B=Point(slices_x[-1], -1, slices_z[-2]),
            C=Point(slices_x[-1], -1, slices_z[-1]),
            D=Point(slices_x[-2], -1, slices_z[-1]),
            ))

vertices = []
vertices_set = set()
vertices_find = {}

faces = []

for square in squares:
    for point in (square.A, square.B, square.C, square.D):
        if point not in vertices_set:
            vertices.append(point)
            vertices_set.add(point)
            vertices_find[point] = len(vertices) - 1

    faces.append((
        vertices_find[square.A],
        vertices_find[square.C],
        vertices_find[square.B]
        ))
    faces.append((
        vertices_find[square.A],
        vertices_find[square.D],
        vertices_find[square.C]
        ))

with open('plane.obj', 'w') as f:
    for vertex in vertices:
        f.write("v %f %f %f\n" % (vertex.x, vertex.y, vertex.z))

    f.write("vt 0 0\n")

    for face in faces:
        f.write("f %d/1 %d/1 %d/1\n" % (face[0] + 1, face[1] + 1, face[2] + 1))

