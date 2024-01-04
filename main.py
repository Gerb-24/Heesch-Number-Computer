from utils import getFirstFromSet

class Vector2:
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y
    
    def equals(self, other):
        return (self.x == other.x and self.y == other.y)
    
    def rotate(self):
        return Vector2(-self.y, self.x)
    
    def flip(self):
        return Vector2(-self.x, self.y)
    
    def add(self, other):
        return Vector2(self.x+other.x, self.y+other.y)
    
    def subtract(self, other):
        return Vector2(self.x-other.x, self.y-other.y)
    
    def __eq__(self, other):
        return isinstance(other, Vector2) and self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))
    
    def __str__(self) -> str:
        return f'Vector2({self.x},{self.y})'

class Tile:
    def __init__(self, coords) -> None:
        self.coords = coords

    def translate(self, translation):
        return Tile(tuple(coord.add(translation) for coord in self.coords))
    
    # computes if two tiles intersect with eachother
    def doIntersect(self, other):
        for coord in self.coords:
            for _coord in other.coords:
                if coord.equals(_coord):
                    return True
        return False
    
    def getBoundary(self, unitVectors):
        counts = []
        boundary = []
        for coord in self.coords:
            for unitVector in unitVectors:
                boundaryCoord = coord.add(unitVector)
                if boundaryCoord in self.coords:
                    continue
                if boundaryCoord in boundary:
                    index = boundary.index(boundaryCoord)
                    counts[index] += 1
                else:
                    boundary.append(boundaryCoord)
                    counts.append(1)
        boundaryWithPriority = createNewBoundary()
        for index, boundaryCoord in enumerate(boundary):
            boundaryWithPriority[counts[index]].add(boundaryCoord)
        return boundaryWithPriority

class CollisionData:
    def __init__(self, collisions=None) -> None:
        self.collisions = collisions if collisions else createNewCollisions()

    def getCollisionsFromBaseOrientations(self, baseTile: "Tile", baseOrientations ):
        for tileType, baseOrientation in baseOrientations.items():
            for orientationCoord in baseOrientation.coords:
                for baseTileCoord in baseTile.coords:
                    translation = baseTileCoord.subtract(orientationCoord)
                    self.collisions[tileType].add( translation )
    
    # we still need to know the collisions for a different shapecode then base shapecode.
    def getFlippedCollisions(self):
        newCollisions = createNewCollisions()
        # when we flip the basetile, then the flipped collisions are for the flipped orientation
        for orientationIndex in self.collisions:
            newCollisions[flipIndex(orientationIndex)] = flipCoords(self.collisions[orientationIndex])
        return CollisionData(newCollisions)
    
    def getRotatedCollisions(self):
        newCollisions = createNewCollisions()
        # when we rotate the basetile, then the rotated collisions are for the rotated orientation
        for orientationIndex in self.collisions:
            newCollisions[rotateIndex(orientationIndex)] = rotateCoords(self.collisions[orientationIndex])
        return CollisionData(newCollisions)

    def getTranslatedCollisions(self, translation):
        newCollisions = createNewCollisions()
        for orientationIndex in self.collisions:
            newCollisions[orientationIndex] = translateCoords(self.collisions[orientationIndex], translation)
        return CollisionData(newCollisions)

    def getCollisionsForBaseOrientations(self):
        collisionsForBaseOrientations = {
        0:      self,
        1:      self.getRotatedCollisions(),
        2:      self.getRotatedCollisions().getRotatedCollisions(),
        3:      self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions(),
        4:      self.getFlippedCollisions(),
        5:      self.getRotatedCollisions().getFlippedCollisions(),
        6:      self.getRotatedCollisions().getRotatedCollisions().getFlippedCollisions(),
        7:      self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions().getFlippedCollisions()
        }
        return collisionsForBaseOrientations


    def getCollisionsFromShapeCode(self, orientationIndex, translation):
        # This works but is not very DRY
        if orientationIndex == 0:
            return self.getTranslatedCollisions(translation)
        elif orientationIndex == 1:
            return self.getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 2:
            return self.getRotatedCollisions().getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 3:
            return self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 4:
            return self.getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 5:
            return self.getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 6:
            return self.getRotatedCollisions().getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)
        elif orientationIndex == 7:
            return self.getRotatedCollisions().getRotatedCollisions().getRotatedCollisions().getFlippedCollisions().getTranslatedCollisions(translation)

class ShapeCode:
    def __init__(self, orientationIndex, translation) -> None:
        self.orientation = orientationIndex
        self.translation = translation

    def __str__(self) -> str:
        return f'Shapecode({self.orientation}, {self.translation})'
    
    def printMotionCanvas(self) -> str:
        return f'[[{self.translation.x},{self.translation.y}],{enumNames[self.orientation]}],'

class Configuration:
    def __init__(self, configuration, boundary) -> None:
        self.configuration = configuration
        self.boundary = boundary

    def getCurrentPriority(self):
        count = 8
        while count > 0:
            if self.boundary[count]:
                return count
            count -= 1
        return 0
    
    def getStartCoord(self):
        startCoord = getFirstFromSet(self.boundary[self.getCurrentPriority()])
        return startCoord

    def getNextConfigurations(self, baseOrientations):
        startCoord = self.getStartCoord()
        validNewTiles =  getValidNewTiles(baseOrientations, self.configuration, startCoord)
        nextConfigurations = set()
        for shapeCode in validNewTiles:
            # If we build up the whole configuration structure as a tree, then the previous configuration is redudant information.
            newConfiguration = Configuration( [*self.configuration, shapeCode ], getNewBoundaryForValidNewTile(self.boundary, baseOrientations, shapeCode) )
            nextConfigurations.add(newConfiguration)
        return nextConfigurations

    def dfs(self, baseOrientations):
        if self.getCurrentPriority() == 0:
            print(self.printMotionCanvas())
            return [self]
        nextConfigurations = self.getNextConfigurations(baseOrientations)
        return [ corona for nextConfiguration in nextConfigurations for corona in nextConfiguration.dfs(baseOrientations) ]
    
    def printMotionCanvas(self):
        baseTileText = self.configuration[0].printMotionCanvas()
        coronaText = '\n'.join([conf.printMotionCanvas() for conf in self.configuration[1:]])
        outText = f'''
  yield* coronas[0].addTiles([
    {baseTileText}
  ], 1)

  const output: [PossibleVector2, TileType][] =  [
    {coronaText}
]
'''
        return outText
        

def rotateCoords( coordsAsVectors ):
    return set(vec.rotate() for vec in coordsAsVectors)

def flipCoords( coordsAsVectors ):
    return set(vec.flip() for vec in coordsAsVectors)

def translateCoords( coordsAsVectors, translate ):
    return set(vec.add(translate) for vec in coordsAsVectors)

def rotateIndex( orientationIndex ):
    if orientationIndex <= 3:
        return (orientationIndex + 1) % 4
    # Apperently we cycle through the rotations the other way when flipped.
    # This should get fixed later?
    return (((orientationIndex - 4) - 1) % 4 ) + 4

def flipIndex( orientationIndex ):
    return (orientationIndex + 4)%8

def createNewCollisions():
    return { index: set() for index in range(amountOfOrientations) }

def createNewBoundary():
    boundaryWithPriority = { i: set() for i in range(1,9)}
    return boundaryWithPriority

amountOfOrientations = 8
allOrientations = ((0,0), (1,0), (2,0), (3,0),
                   (0,1), (1,1), (2,1), (3,1))
enumNames = ('TileType.F0', 'TileType.F90', 'TileType.F180', 'TileType.F270',
             'TileType.T0', 'TileType.T90', 'TileType.T180', 'TileType.T270')

unitVectors = (Vector2(1,0), Vector2(-1,0), Vector2(0,1), Vector2(0, -1), Vector2(1,1), Vector2(-1,-1), Vector2(1,-1), Vector2(-1,1))
# coords = [[0, 0],
# 		[1, 0], 
# 		[2, 0], [2, 1], [2, 2],]
coords =  [[0, 0],
		[1, 0], [1, -1],
		[2, -1], [2, -3],
		[3, 0], [3, -1], [3, -2], [3, -3],
		[4, -3],
		[5, -2], [5, -3],
		[6, -3]]


coordsAsVectors = set(Vector2(*coord) for coord in coords)


# Here we construct the tiles with the different rotations and amount of flips.
baseOrientations = {
    0:      Tile(coordsAsVectors),
    1:      Tile(rotateCoords(coordsAsVectors)),
    2:      Tile(rotateCoords(rotateCoords(coordsAsVectors))),
    3:      Tile(rotateCoords(rotateCoords(rotateCoords(coordsAsVectors)))),
    4:      Tile(flipCoords(coordsAsVectors)),
    5:      Tile(flipCoords(rotateCoords(coordsAsVectors))),
    6:      Tile(flipCoords(rotateCoords(rotateCoords(coordsAsVectors)))),
    7:      Tile(flipCoords(rotateCoords(rotateCoords(rotateCoords(coordsAsVectors)))))
}

baseTile = baseOrientations[0]
baseBoundary = baseTile.getBoundary(unitVectors)

# using these we compute the different translations in which these intersect with out base tile. 
# We can then forget about the tiles coordinates, as this is the only data that matters.

collisionData = CollisionData()
collisionData.getCollisionsFromBaseOrientations( baseTile, baseOrientations )
collisionsForBaseOrientations = collisionData.getCollisionsForBaseOrientations()

print('all collisions are computed')

# For the first corona we take loop over the coordinates in the boundary and want to fill them in.
# So we want to loop over all the coordinates of all the baseOrientations, 
# and see if the translation to the corresponding square in the boundary 
# is not in the set of collisions for all of the tiles in the current configuration set. 
# If it is not, then we make a deepcopy of our configuration set, where we add this new tile.

# The baseConfiguration is a list of ShapeCodes,
# we then want to find all the new ShapeCodes so that they cover the startCoord but not the  tiles in the baseConfiguration

def getValidNewTiles(baseOrientations, configuration, startCoord: Vector2):
    validNewTiles = set()
    for orientationIndex, baseOrientation in baseOrientations.items():
        for coord in baseOrientation.coords:
            doesNotCollide = True
            for shapeCodeInConfiguration in configuration:
                collisionsForShapeCodeInConfiguration = collisionsForBaseOrientations[shapeCodeInConfiguration.orientation].getTranslatedCollisions(shapeCodeInConfiguration.translation)
                collisionsFromShapeCodeForOrientationIndex = collisionsForShapeCodeInConfiguration.collisions[orientationIndex]
                translation = startCoord.subtract(coord)
                if translation in collisionsFromShapeCodeForOrientationIndex:
                    doesNotCollide = False
                    break
            if doesNotCollide:
                validNewTiles.add( ShapeCode(orientationIndex, translation) )
    return validNewTiles

# We now have the valid new tiles, so we add these to get a new configuration.
# To go on with adding tiles, we need to find the new startingCoord.
# Going to the next tile in the baseBoundary does not work, our added tile might already cover it.
# So we need to update our boundary and continue to use this

def getTilesFromShapeCode(baseOrientations, shapeCode: ShapeCode):
    baseOrientation = baseOrientations[ shapeCode.orientation ]
    tiles = translateCoords(baseOrientation.coords, shapeCode.translation)
    return tiles

def getNewBoundaryForValidNewTile(currentBoundary, baseOrientations, shapeCode: ShapeCode):
    newTiles = getTilesFromShapeCode(baseOrientations, shapeCode)
    newBoundary = createNewBoundary()
    for count in currentBoundary:
        newBoundary[count] = currentBoundary[count].difference(newTiles)
    return newBoundary


baseConfiguration = Configuration( [ShapeCode(0, Vector2(0,0))], baseBoundary)
print('starting now!')
ans = baseConfiguration.dfs(baseOrientations)
print('finished!')