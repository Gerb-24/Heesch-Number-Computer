from squareshapes import Square, Polyomino
from main import Vector2

coords = ([0,0],
     [1,0], [1,-1], [1,1], 
     [2,1], 
     [3,1], [3,2], 
     [4,0], [4,1], [4,2], 
     [5,0])

coords2 = [[0, 0],
		[1, 0], [1, -1],
		[2, -1], [2, -3],
		[3, 0], [3, -1], [3, -2], [3, -3],
		[4, -3],
		[5, -2], [5, -3],
		[6, -3]]

coords =  [[0, 0],
		[1, 0],
		[2, 0], [2, 1],[2,2],
		[3, 1],
		[4, 1], [4,0], [4,-1], [4,-2],
		[5, -1],
		[6, -1], [6,0]]

coords = [[0, 0], [0, 1], [0, 2],
		[1, 0], [1, 2],
		[2, 0], [2, 1], [2, 2],]

coords =     [[0, 0],
		[1, 0], [1,1], [1,2], [1,3],
         [2,2] ]

# 2.5 secs
tile = Polyomino(
	[
		Square(0, 0),
		Square(2, 0), Square(2, 2),
		Square(4, 2), Square(4, 6),
		Square(6, 0), Square(6, 2), Square(6, 4), Square(6, 6),
		Square(8, 6),
		Square(10, 4), Square(10, 6),
		Square(12, 6),
	],
	priority = [
		Square(4,4),
		Square(8,4),
		Square(12,4),
		Square(8,2),
		Square(4,0),
	]
)

# 2.5 secs
coords =  [[0, 0],
		[1, 0], [1, -1],
		[2, -1], [2, -3],
		[3, 0], [3, -1], [3, -2], [3, -3],
		[4, -3],
		[5, -2], [5, -3],
		[6, -3]]
priority = [
		Vector2(2,-2),
		Vector2(4,-2),
		Vector2(6,-2),
		Vector2(4,-1),
		Vector2(2,0),
	]

# 22 secs
tile = Polyomino(
	[
		Square(0, 4),
		Square(2, 0), Square(2, 2), Square(2, 4),
		Square(4, 0), Square(4, 2), Square(4, 4),
		Square(6, 4),
		Square(8, 4), Square(8, 6),
		Square(10, 6),
	],
	priority = [
		Square(6,6),
		Square(10,4),
		Square(0,2),
		Square(6,2),
	]
)

# 22 secs
coords = [
		(0, -2),
		(1, 0), (1, -1), (1, -2),
		(2, 0), (2, -1), (2, -2),
		(3, -2),
		(4, -2), (4, -3),
		(5, -3),
	]
priority = [
	Vector2(3,-3),
	Vector2(5,-2),
	Vector2(0,-1),
	Vector2(3,-1),
]


# 90 secs
tile = Polyomino(
	[
		Square(0, 4),
		Square(2, 4), Square(2, 6), Square(2, 8),
		Square(4, 0), Square(4, 2), Square(4, 4),
		Square(6, 0), Square(6, 2),
		Square(8, 0), Square(8, 2),
	],
	priority = [
		Square(0,6),
		Square(4,6),
		Square(6,4),
		Square(2,2),
	]
)