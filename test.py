from squareshapes import Square, Polyomino
import time

start_time = time.time()
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

heesch_configs = tile.heesch_computer()
print("--- %s seconds ---" % (time.time() - start_time))