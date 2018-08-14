using SixLabors.ImageSharp;
using SixLabors.ImageSharp.PixelFormats;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace Planetary_Generation
{
    /// <summary>
    /// Contains data and methods for managing the entire map.
    /// </summary>
    public class PlateLayer
    {
        /// <summary>
        /// Number of points per height of map, twice this number for length of map.
        /// </summary>
        private static int _layerSize;
        /// <summary>
        /// Contains all the plates on the map, independently from <see cref="Map"/>.
        /// </summary>
        private readonly Plate[] Plates;
        /// <summary>
        /// Contains all the points on the map, independently from <see cref="Plates"/>.
        /// </summary>
        private readonly Point[,] Map;
        /// <summary>
        /// Contains extra lists for keeping track of points during some methods.
        /// </summary>
        private readonly Plate testPlate;
        /// <summary>
        /// Matrix of points including which plate each point is part of.
        /// </summary>
        private readonly int[,] plateTracker;
        /// <summary>
        /// Allocates space for arrays and sets up the point class.
        /// </summary>
        /// <param name="iSize">Size of the map, see <see cref="_layerSize"/></param>
        /// <param name="iNumberPlates">Number of plates in <see cref="Plates"/></param>
        public PlateLayer(int iSize, int iNumberPlates)
        {
            _layerSize = iSize;
            Plates = new Plate[iNumberPlates];
            for (int i = 0; i < Plates.Length; i++)
            {
                Plates[i] = new Plate();
            }
            Map = new Point[2 * _layerSize, _layerSize];
            Point.MapSetup(iSize);
            for (int x = 0; x < 2 * _layerSize; x++)
            {
                for (int y = 0; y < _layerSize; y++)
                {
                    Map[x, y] = new Point(x, y);
                }
            }
            testPlate = new Plate();
            plateTracker = new int[2 * _layerSize, _layerSize];
        }
        /// <summary>
        /// Increases the Height of all points within circles centered at points with <see cref="Point.IsMember"/> set to true.
        /// </summary>
        /// <param name="radius">Radius of circles.</param>
        /// <param name="magnitude">Magnitude of height to be added per point per circle.</param>
        public void FlowStep(double radius, double magnitude)
        {
            double radiusSquared = radius * radius;
            Parallel.For(0, (2 * _layerSize), (i) =>
            {
                for (int j = 0; j < _layerSize; j++)
                {
                    if (Map[i, j].IsMember)
                    {
                        Map[i, j].IsMember = false;
                        int[] keyRange = new int[4];
                        keyRange = Map[i, j].Range(radius);
                        for (int xi = keyRange[0]; xi < keyRange[1]; xi++)
                        {
                            for (int y = keyRange[2]; y < keyRange[3]; y++)
                            {
                                int x = xi;
                                if (x >= 2 * _layerSize)
                                {
                                    x -= 2 * _layerSize;
                                }
                                if (Map[i, j].Distance(Map[x, y]) < radiusSquared)
                                {
                                    Map[x, y].Height += magnitude;
                                }
                            }
                        }
                    }
                }
            });
        }
        /// <summary>
        /// Sorts all points and returns the height at the given index.
        /// </summary>
        /// <param name="index">Index to return height value.</param>
        /// <returns>Height at given index.</returns>
        public double MediumHeight(int index)
        {
            double[] unsorted = new double[2 * _layerSize * _layerSize];
            int k = 0;
            foreach (Point iPoint in Map)
            {
                unsorted[k] = iPoint.Height;
                k++;
            }
            Array.Sort(unsorted);
            return unsorted[index];
        }
        /// <summary>
        /// Sets <see cref="Point.IsMember"/> to true for points above height threshold, false otherwise.
        /// </summary>
        /// <param name="midHeight">Threshold height.</param>
        public void HeightFilterAlpha(double midHeight)
        {
            foreach (Point iPoint in Map)
            {
                if (midHeight < iPoint.Height)
                {
                    iPoint.IsMember = true;
                }
                else
                {
                    iPoint.IsMember = false;
                }
                iPoint.Height = 0;
            }
        }
        /// <summary>
        /// Scans for all points in <see cref="Map"/> contiguously adjacent to the given point with <see cref="Point.IsMember"/>
        /// set to true and adds them to <see cref="testPlate"/>, setting <see cref="Point.IsMember"/> to false in the process.
        /// </summary>
        /// <param name="x">XCoordinate of given point.</param>
        /// <param name="y">YCoordinate of given point.</param>
        public void CheckNeighbor(int x, int y)
        {
            Stack<int[]> pointStack = new Stack<int[]>();
            int[] pointA = { x, y };
            Map[x, y].IsMember = false;
            pointStack.Push(pointA);
            while (pointStack.Count != 0)
            {
                int[] point = pointStack.Pop();
                testPlate.Value.Add(new Point(point[0], point[1]));
                for (int i = 0; i < 4; i++)
                {
                    int[] pointB = new int[2];
                    pointB[0] = Map[point[0], point[1]].NeighborPoints[0, i];
                    pointB[1] = Map[point[0], point[1]].NeighborPoints[1, i];
                    if (Map[pointB[0], pointB[1]].IsMember)
                    {
                        Map[pointB[0], pointB[1]].IsMember = false;
                        pointStack.Push(pointB);
                    }
                }
            }
        }
        /// <summary>
        /// Uses <see cref="CheckNeighbor(int, int)"/> to generate a plate starting at the given point, and transfers it
        /// to <see cref="Plates"/> if there is an empty plate or the plate is larger than the smallest existing plate.
        /// </summary>
        /// <param name="x">XCoordinate of given point.</param>
        /// <param name="y">YCoordinate of given point.</param>
        public void PlateMaker(int x, int y)
        {
            if (testPlate.Value != null)
            {
                testPlate.Value.Clear();
            }
            if (Map[x,y].IsMember)
            {
                CheckNeighbor(x, y);
            }
            foreach (Plate iPlate in Plates)
            {
                if (iPlate.Value.Count == 0)
                {
                    iPlate.Copy(testPlate);
                    return;
                }
            }
            foreach (Plate iPlate in Plates)
            {
                if (iPlate.Value.Count < testPlate.Value.Count)
                {
                    iPlate.Copy(testPlate);
                    return;
                }
            }
        }
        /// <summary>
        /// Generates the largest <see cref="Plates"/>, sets <see cref="Point.IsMember"/> to true for every point in one
        /// of those plates, and false otherwise.
        /// </summary>
        public void PlateMaking()
        {
            for (int x = 0; x < 2 * _layerSize; x++)
            {
                for (int y = 0; y < _layerSize; y++)
                {
                    if (Map[x, y].IsMember)
                    {
                        PlateMaker(x, y);
                    }
                }
            }
            foreach (Plate iPlate in Plates)
            {
                foreach (Point iPoint in iPlate.Value)
                {
                    Map[iPoint.XCoord, iPoint.YCoord].IsMember = true;
                }
            }
        }
        /// <summary>
        /// Expands each plate one pixel in every direction until no point outside of all plates exists.
        /// </summary>
        public void ExpandPlates()
        {
            Queue<int[]> borderPoints = new Queue<int[]>();
            for (int i = 0; i < Plates.Length; i++)
            {
                for (int j = 0; j < Plates[i].Value.Count; j++)
                {
                    int xCoord = Plates[i].Value[j].XCoord;
                    int yCoord = Plates[i].Value[j].YCoord;
                    for (int k = 0; k < 4; k++)
                    {
                        int xCoordA = Map[xCoord, yCoord].NeighborPoints[0, k];
                        int yCoordA = Map[xCoord, yCoord].NeighborPoints[1, k];
                        if (!Map[xCoordA, yCoordA].IsMember)
                        {
                            int[] borderPoint = new int[3] {xCoordA, yCoordA, i};
                            borderPoints.Enqueue(borderPoint);
                        }
                    }
                }
            }
            while (borderPoints.Count != 0)
            {
                int[] index = borderPoints.Dequeue();
                if (!Map[index[0],index[1]].IsMember)
                {
                    Plates[index[2]].Value.Add(new Point(index[0], index[1]));
                    Map[index[0], index[1]].IsMember = true;
                    for (int k = 0; k < 4; k++)
                    {
                        int xCoordA = Map[index[0], index[1]].NeighborPoints[0, k];
                        int yCoordA = Map[index[0], index[1]].NeighborPoints[1, k];
                        if (!Map[xCoordA, yCoordA].IsMember)
                        {
                            int[] borderPoint = new int[3] { xCoordA, yCoordA, index[2]};
                            borderPoints.Enqueue(borderPoint);
                        }
                    }
                }
            }
        }
        /// <summary>
        /// Generates plates based off of sets of circles with given probabilistic distribution corrected for area. 
        /// Then, expands plates to include all points.
        /// </summary>
        /// <param name="magnitude">Array of weight of each set of circles.</param>
        /// <param name="radius">Radius of each set of circles.</param>
        /// <param name="pointConcentration">Raw probability that a point will be in a circle.</param>
        /// <param name="cutOff">Approximate number of points which will initially become plates.</param>
        public void PlateGeneration(double[] magnitude, double[] radius, double[] pointConcentration, int cutOff)
        {
            Random rnd = new Random();
            for (int i = 0; i < magnitude.Length; i++)
            {
                foreach (Point iPoint in Map)
                {
                    iPoint.RandomMomentum(rnd.NextDouble(), pointConcentration[i]);
                }
                FlowStep(radius[i], magnitude[i]);
            }
            HeightFilterAlpha(MediumHeight(cutOff));
            PlateMaking();
            ExpandPlates();
        }
        /// <summary>
        /// Saves data containing which plate each point belongs to, to a BMP image file.
        /// </summary>
        /// <param name="directory">Folder to contain the file. Does not include final /.</param>
        /// <param name="fileName">Name of file. Does not include .bmp extension.</param>
        public void SavePlates(string directory, string fileName)
        {
            double sd = (510 / (Plates.Length - 1));
            int si;
            using (Image<Rgba32> image = new Image<Rgba32>(2 * _layerSize, _layerSize))
            {
                for (int n = 0; n < Plates.Length; n++)
                {
                    Rgba32 plateColor;
                    plateColor.A = 255;
                    si = (int)Math.Round((double)n * sd) - 255;
                    if (si < 0)
                    {
                        plateColor.R = (byte)(-1 * si);
                        plateColor.G = 0;
                        plateColor.B = 0;
                    }
                    else
                    {
                        plateColor.R = 0;
                        plateColor.G = (byte)si;
                        plateColor.B = (byte)si;
                    }
                    foreach (Point iPoint in Plates[n].Value)
                    {
                        image[iPoint.XCoord, iPoint.YCoord] = plateColor;
                    }
                }
                image.Save(directory + "\\" + fileName + ".png");
            }
        }
        /// <summary>
        /// Saves the height of each point to a BMP image file.
        /// </summary>
        /// <param name="directory">Folder to contain the file. Does not include final /.</param>
        /// <param name="fileName">Name of file. Does not include .bmp extension.</param>
        public void SaveHeights(string directory, string fileName)
        {
            double minHeight = double.MaxValue;
            double maxHeight = 0;
            foreach (Point iPoint in Map)
            {
                if (iPoint.Height < minHeight)
                {
                    minHeight = iPoint.Height;
                }
                if (iPoint.Height > maxHeight)
                {
                    maxHeight = iPoint.Height;
                }
            }
            maxHeight = (510 / (maxHeight - minHeight));
            using (Image<Rgba32> image = new Image<Rgba32>(2 * _layerSize, _layerSize))
            {
                foreach (Point iPoint in Map)
                {
                    Rgba32 plateColor;
                    plateColor.A = 255;
                    int si = (int)Math.Round((iPoint.Height - minHeight) * maxHeight) - 255;
                    if (si < 0)
                    {
                        plateColor.R = (byte)(-1*si);
                        plateColor.G = 0;
                        plateColor.B = 0;
                    }
                    else
                    {
                        plateColor.R = 0;
                        plateColor.G = (byte)si;
                        plateColor.B = (byte)si;
                    }
                    image[iPoint.XCoord, iPoint.YCoord] = plateColor;
                }
                image.Save(directory + "\\" + fileName + ".png");
            }
        }
        /// <summary>
        /// Stores which points belonged to which plate in <see cref="plateTracker"/>.
        /// </summary>
        private void PlateTracker()
        {
            for (int i = 0; i < Plates.Length; i++)
            {
                foreach (Point iPoint in Plates[i].Value)
                {
                    plateTracker[iPoint.XCoord, iPoint.YCoord] = i;
                }
            }
        }
        /// <summary>
        /// Calculates the average height of adjacent points in <see cref="Map"/>.
        /// </summary>
        /// <param name="iPlate">Index of plate to scan for with <see cref="plateTracker"/>.</param>
        /// <param name="x">X Coordinate of center point.</param>
        /// <param name="y">Y Coordinate of center point.</param>
        /// <returns>Average of nearby points from <see cref="Map"/>.</returns>
        private double AdjacentHeightAverage(int iPlate, int x, int y, bool moveHorizontal)
        {
            Point tempPoint = new Point(x, y);
            int[,] neighborPoints = tempPoint.NeighborPoints;
            double average = 0;
            int averageCount = 0;
            if (moveHorizontal)
            {
                for (int i = 2; i < 4; i++)
                {
                    if (plateTracker[neighborPoints[0, i], neighborPoints[1, i]] == iPlate)
                    {
                        if (Map[neighborPoints[0, i], neighborPoints[1, i]].Height != 0)
                        {
                            average += Map[neighborPoints[0, i], neighborPoints[1, i]].Height;
                            averageCount++;
                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < 2; i++)
                {
                    if (plateTracker[neighborPoints[0, i], neighborPoints[1, i]] == iPlate)
                    {
                        if (Map[neighborPoints[0, i], neighborPoints[1, i]].Height != 0)
                        {
                            average += Map[neighborPoints[0, i], neighborPoints[1, i]].Height;
                            averageCount++;
                        }
                    }
                }
            }
            if (averageCount == 0)
            {
                for (int i = 0; i < 4; i++)
                {
                    if (plateTracker[neighborPoints[0, i], neighborPoints[1, i]] == iPlate)
                    {
                        if (Map[neighborPoints[0, i], neighborPoints[1, i]].Height != 0)
                        {
                            average += Map[neighborPoints[0, i], neighborPoints[1, i]].Height;
                            averageCount++;
                        }
                    }
                }
            }
            if (averageCount > 1)
            {
                average = average / (double)averageCount;
            }
            return average;
        }
        /// <summary>
        /// Adds points to plates by tracing them backwards.
        /// </summary>
        /// <param name="input">Point to add.</param>
        /// <param name="timeStep">Time factor for moving plates.</param>
        private void ReverseAdd(Point input, double timeStep)
        {
            for (int i = 0; i < Plates.Length; i++)
            {
                int[] reversedCoord = new int[2];
                double[] angle = new double[3]
                {
                    Plates[i].Direction[0], Plates[i].Direction[1], -1 * Plates[i].Speed * timeStep
                };
                reversedCoord = input.Transform(angle);
                if (plateTracker[reversedCoord[0], reversedCoord[1]] == i)
                {
                    bool horizontalMove = false;
                    if (input.YCoord == reversedCoord[1])
                    {
                        horizontalMove = true;
                    }
                    double average = AdjacentHeightAverage(i, reversedCoord[0], reversedCoord[1], horizontalMove);
                    Plates[i].Value.Add(new Point(input.XCoord, input.YCoord, average));
                    Map[input.XCoord, input.YCoord].IsMember = true;
                }
            }
        }
        /// <summary>
        /// Small Class for sorting arrays of integers by the first two entries.
        /// </summary>
        public class FourPointComp : IComparer<int[]>
        {
            public int Compare(int[] x, int[] y)
            {
                if (x is int[] xi)
                {
                    if (y is int[] yi)
                    {
                        if (xi[1].CompareTo(yi[1]) != 0)
                        {
                            return xi[1].CompareTo(yi[1]);
                        }
                        if (xi[0].CompareTo(yi[0]) != 0)
                        {
                            return xi[0].CompareTo(yi[0]);
                        }
                    }
                }
                return 0;
            }
        }
        /// <summary>
        /// Finds all points with two or more entries in multiple plates and adds them to <see cref="testPlate"/>.
        /// </summary>
        private List<List<int[]>> FindOverlap()
        {
            SortedSet<int[]> output = new SortedSet<int[]>(new FourPointComp());
            bool[,] MapTwo = new bool[2 * _layerSize, _layerSize];
            Parallel.For(0, Plates.Length, (i) =>
            {
                foreach (Point iPoint in Plates[i].Value)
                {
                    if (Map[iPoint.XCoord, iPoint.YCoord].IsMember)
                    {
                        Map[iPoint.XCoord, iPoint.YCoord].IsMember = false;
                    }
                    else
                    {
                        MapTwo[iPoint.XCoord, iPoint.YCoord] = true;
                    }
                }
            });
            for (int i = 0; i < Plates.Length; i++)
            {
                for (int j = 0; j < Plates[i].Value.Count; j++)
                {
                    if (MapTwo[Plates[i].Value[j].XCoord, Plates[i].Value[j].YCoord])
                    {
                        MapTwo[Plates[i].Value[j].XCoord, Plates[i].Value[j].YCoord] = true;
                        int[] outPoint = new int[4]
                        {Plates[i].Value[j].XCoord, Plates[i].Value[j].YCoord, i, j };
                        output.Add(outPoint);
                    }
                }
            }
            int index = 0;
            int xPt = output.First()[0];
            int yPt = output.First()[1];
            List<List<int[]>> outputList = new List<List<int[]>>();
            outputList.Add(new List<int[]>());
            while (output.Count != 0)
            {
                if (xPt == output.First()[0] && yPt == output.First()[1])
                {
                    outputList[index].Add(output.First());
                    output.Remove(output.First());
                }
                else
                {
                    index++;
                    outputList.Add(new List<int[]>());
                    xPt = output.First()[0];
                    yPt = output.First()[1];
                }
            }
            return outputList;
        }
        /// <summary>
        /// Corrects to one plate entry for a given point.
        /// </summary>
        /// <param name="input">Input points.</param>
        private void ResolveOverlap(List<int[]> input)
        {
            int[] point = input.First();
            int plateIndex = input[0][2];
            int pointIndex = input[0][3];
            for (int i = 1; i < input.Count; i++)
            {
                if (Plates[input[i][2]].Value[input[i][3]].Height != 0)
                {
                    if (Plates[plateIndex].Value[pointIndex].Height != 0)
                    {
                        Plates[plateIndex].Value[pointIndex].Height += 0.6 * Plates[input[i][2]].Value[input[i][3]].Height;
                        Plates[input[i][2]].Value.RemoveAt(input[i][3]);
                    }
                    else
                    {
                        Plates[plateIndex].Value.RemoveAt(pointIndex);
                        plateIndex = input[i][2];
                        pointIndex = input[i][3];
                    }
                }
                else
                {
                    Plates[input[i][2]].Value.RemoveAt(input[i][1]);
                }
            }
        }
        /// <summary>
        /// Moves all plates for a given time factor.
        /// </summary>
        /// <param name="timeStep">Time factor for plate movement.</param>
        public void PlateMove(double timeStep)
        {
            PlateTracker();
            Parallel.For(0, Plates.Length, (i) =>
            {
                Plates[i].Slide(timeStep);
            });
            Parallel.For(0, Plates.Length, (i) =>
            {
                foreach (Point iPoint in Plates[i].Value)
                {
                    Map[iPoint.XCoord, iPoint.YCoord].IsMember = true;
                }
            });
            foreach (Point iPoint in Map)
            {
                if (!iPoint.IsMember)
                {
                    ReverseAdd(iPoint, timeStep);
                }
            }
            ExpandPlates();
            Parallel.For(0, Plates.Length, (i) =>
            {
                Plates[i].CondensePoints();
            });
            List<List<int[]>> overlapList =  FindOverlap();
            foreach (List<int[]> iList in overlapList)
            {
                ResolveOverlap(iList);
            }
            Parallel.For(0, Plates.Length, (i) =>
            {
                foreach (Point iPoint in Plates[i].Value)
                {
                    Map[iPoint.XCoord, iPoint.YCoord].Height = iPoint.Height;
                }
            });
        }
        /// <summary>
        /// Inputs height data from a given BMP image file.
        /// </summary>
        /// <param name="directory">Folder containing image file, not including final /.</param>
        /// <param name="heightFile">Image file name, not including .bmp extension.</param>
        /// <returns>True if successful, false otherwise.</returns>
        public bool HeightInputs(string directory, string heightFile)
        {
            Rgba32 pixelColor;
            try
            {
                using (Image<Rgba32> image = Image.Load(directory + "\\" + heightFile + ".png"))
                {
                    for (int x = 0; x < 2 * _layerSize; x++)
                    {
                        for (int y = 0; y < _layerSize; y++)
                        {
                            pixelColor = image[x, y];
                            Map[x, y].Height = (pixelColor.G - 1*pixelColor.R) + 255;
                        }
                    }
                }
            }
            catch (Exception)
            {
                return false;
            }
            return true;
        }
        /// <summary>
        /// Inputs which points belong to which plates from a given BMP image file.
        /// </summary>
        /// <param name="directory">Folder containing image file, not including final /.</param>
        /// <param name="heightFile">Image file name, not including .bmp extension.</param>
        /// <returns>True if successful, false otherwise.</returns>
        public bool PlatesMapInputs(string directory, string plateFile)
        {
            Rgba32 pixelColor;
            try
            {
                List<int> plateIndex = new List<int>();
                int testColor;
                bool addedPlate = false;
                using (Image<Rgba32> image = Image.Load(directory + "\\" + plateFile + ".png"))
                {
                    for (int x = 0; x < 2 * _layerSize; x++)
                    {
                        for (int y = 0; y < _layerSize; y++)
                        {
                            addedPlate = false;
                            pixelColor = image[x, y];
                            testColor = (pixelColor.G - 1 * pixelColor.R) + 255;
                            foreach (int iInt in plateIndex)
                            {
                                if (testColor == iInt)
                                {
                                    addedPlate = true;
                                    break;
                                }
                            }
                            if (!addedPlate)
                            {
                                plateIndex.Add(testColor);
                            }
                        }
                    }
                    if (plateIndex.Count != Plates.Length)
                    {
                        return false;
                    }
                    plateIndex.Sort();
                    for (int x = 0; x < 2 * _layerSize; x++)
                    {
                        for (int y = 0; y < _layerSize; y++)
                        {
                            pixelColor = image[x, y];
                            for (int i = 0; i < Plates.Length; i++)
                            {
                                if (((pixelColor.G - 1 * pixelColor.R) + 255) == plateIndex[i])
                                {
                                    Plates[i].Value.Add(new Point(x, y, Map[x, y].Height));
                                    break;
                                }
                            }
                        }
                    }
                }
                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }
        /// <summary>
        /// Inputs speed and direction data.
        /// </summary>
        /// <param name="speeds">Speed of plate movement to input to <see cref="Plate.Speed"/>.</param>
        /// <param name="directions">Direction of plate movement to input to <see cref="Plate.Direction"/>.</param>
        public void PlatesVelocityInputs(double[] speeds, double[,] directions)
        {
            for (int i = 0; i < Plates.Length; i++)
            {
                Plates[i].Speed = speeds[i];
                Plates[i].Direction[0] = directions[0, i];
                Plates[i].Direction[1] = directions[1, i];
            }
        }
    }
}