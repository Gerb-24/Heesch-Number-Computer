from squareshapes import Square, Polyomino, test_check_combinations
import time

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

start_time = time.time()
# ans = tile.corona_maker(tile.orientations())
tile.heesch_computer()
print("--- %s seconds ---" % (time.time() - start_time))

# for square in tile.outside():
#     print(square)

# for i in range(200):
# 	test_check_combinations(tile, tile.orientations(), [tile], Square(8,4))
# print("--- %s seconds ---" % (time.time() - start_time))
# for tile in test_check_combinations(tile, tile.orientations(), [tile], Square(8,4)):
# 	print(tile.shapecode)