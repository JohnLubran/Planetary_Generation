using System;

namespace Planetary_Generation
{
    /// <summary>
    /// Contains data and functions for points.
    /// </summary>
    public class Point : IComparable
    {
        /// <summary>
        /// Number of points per height of map, twice this number for length of map.
        /// </summary>
        private static int _mapSize;
        /// <summary>
        /// Angular width and height of point in radians.
        /// </summary>
        private static double _angleSize;
        /// <summary>
        /// Constant for converting between array index (YCoord) and phi.
        /// </summary>
        private static double _phiShift;
        /// <summary>
        /// X value of point, for array index.
        /// </summary>
        public readonly int XCoord;
        /// <summary>
        /// Y value of point, for array index.
        /// </summary>
        public readonly int YCoord;
        /// <summary>
        /// Contains index (XCoord,YCoord) of adjacent points.
        /// <para>
        /// ([0,0],[1,0]) and ([0,1],[1,1]) are right and left.
        /// ([0,2],[1,2]) and ([0,3],[1,3]) are above and below.
        /// </para>
        /// </summary>
        public int[,] NeighborPoints;
        /// <summary>
        /// Bool used to keep track of which points have been processed for a number of methods.
        /// </summary>
        public bool IsMember;
        /// <summary>
        /// Height at that point, both for plate movement and plate generation.
        /// </summary>
        public double Height;
        /// <summary>
        /// Angular position of point, ranges from 0 to 2 pi.
        /// </summary>
        private readonly double _theta;
        /// <summary>
        /// Cosine of angular position of point. Phi ranges from -pi/2 to pi/2.
        /// </summary>
        private readonly double _cosPhi;
        /// <summary>
        /// Sine of angular position of point. Phi ranges from -pi/2 to pi/2.
        /// </summary>
        private readonly double _sinPhi;
        /// <summary>
        /// Sets the static variables for the class.
        /// </summary>
        /// <param name="size">Size of map. Will become <see cref="_mapSize"/></param>
        public static void MapSetup(int size)
        {
            _mapSize = size;
            _angleSize = Math.PI / (double)_mapSize;
            _phiShift = 0.5 * _angleSize - (Math.PI / 2);
        }
        /// <summary>
        /// Calculates and stores data for all internal fields.
        /// </summary>
        /// <param name="xIn">XCoord of point, see <see cref="XCoord"/></param>
        /// <param name="yIn">YCoord of point, see <see cref="YCoord"/></param>
        /// <param name="heightIn">Height of point, see <see cref="Height"/></param>
        public Point(int xIn, int yIn, double heightIn = 0)
        {
            XCoord = xIn;
            YCoord = yIn;
            Height = heightIn;
            _theta = XCoord * _angleSize;
            _cosPhi = Math.Cos((YCoord * _angleSize) + _phiShift);
            _sinPhi = Math.Sin((YCoord * _angleSize) + _phiShift);
            IsMember = false;
            NeighborPoints = new int[2, 4];
            FindNeighborPoints(XCoord, YCoord);
        }
        /// <summary>
        /// Determines adjacent points, including wrap-around on the horizontal and vertical boundaries.
        /// </summary>
        /// <param name="x">XCoord of point, see <see cref="XCoord"/></param>
        /// <param name="y">YCoord of point, see <see cref="YCoord"/></param>
        private void FindNeighborPoints(int x, int y)
        {
            NeighborPoints[1, 0] = y;
            NeighborPoints[1, 1] = y;
            if (x == 0)
            {
                NeighborPoints[0, 0] = 2 * _mapSize - 1;
                NeighborPoints[0, 1] = 1;
            }
            else if (x == 2 * _mapSize - 1)
            {
                NeighborPoints[0, 0] = 0;
                NeighborPoints[0, 1] = 2 * _mapSize - 2;
            }
            else
            {
                NeighborPoints[0, 0] = x + 1;
                NeighborPoints[0, 1] = x - 1;
            }
            NeighborPoints[0, 2] = x;
            if (y == 0 || y == _mapSize - 1)
            {
                int newI = x + _mapSize;
                if (newI >= 2 * _mapSize)
                {
                    newI = newI - 2 * _mapSize;
                }
                if (y == 0)
                {
                    NeighborPoints[1, 2] = y + 1;
                }
                else
                {
                    NeighborPoints[1, 2] = y - 1;
                }
                NeighborPoints[0, 3] = newI;
                NeighborPoints[1, 3] = y;
            }
            else
            {
                NeighborPoints[1, 2] = y + 1;
                NeighborPoints[0, 3] = x;
                NeighborPoints[1, 3] = y - 1;
            }
        }
        /// <summary>
        /// Scales the random probability for area, and sets <see cref="IsMember"/> for points that exceed the threshold.
        /// </summary>
        /// <param name="randomNumber">Random number, unscaled for point size.</param>
        /// <param name="threshhold">Minimum value to pass probability test.</param>
        public void RandomMomentum(double randomNumber, double threshhold)
        {
            if (randomNumber > Math.Pow(threshhold, _cosPhi))
            {
                IsMember = true;
            }
        }
        /// <summary>
        /// Calculates rectangular approximation of points in range of this point.
        /// </summary>
        /// <param name="range">Distance from center point.</param>
        /// <returns>Limits of points in range.
        /// [0] to [1] are X limits, adding 2*<see cref="_mapSize"/> if range wraps around border.
        /// [2] to [3] are y limits.</returns>
        public int[] Range(double range)
        {
            int[] output = new int[4];
            double zMin = ((_sinPhi - range) + 1) / 2;
            if (zMin < 0)
            {
                zMin = 0;
            }
            double zMax = ((_sinPhi + range) + 1) / 2;
            if (zMax > 1)
            {
                zMax = 1;
            }
            double rDif = range / _cosPhi;
            if (rDif > 1)
            {
                rDif = 1;
            }
            output[0] = XCoord - (int)(rDif * (double)_mapSize);
            output[1] = XCoord + (int)(rDif * (double)_mapSize);
            if (output[0] <= 0)
            {
                output[0] += 2 * _mapSize;
                output[1] += 2 * _mapSize;
            }
            output[2] = (int)Math.Floor(zMin * (double)_mapSize);
            output[3] = (int)Math.Ceiling(zMax * (double)_mapSize);
            return output;
        }
        /// <summary>
        /// Calculates distance between this point and the input point.
        /// </summary>
        /// <param name="iPoint">Second point used for calculating distance.</param>
        /// <returns>Distance between the two points.</returns>
        public double Distance(Point iPoint)
        {
            int xDif = Math.Abs(iPoint.XCoord - XCoord);
            if (xDif > _mapSize)
            {
                xDif = 2 * _mapSize - xDif;
            }
            double cosDif;
            if (xDif < _mapSize/6)
            {
                cosDif = 1 - (0.5 * xDif * xDif* _angleSize * _angleSize);
            }
            else
            {
                cosDif = Math.Cos((double)xDif * _angleSize);
            }
            return Math.Abs(2*(1 - _sinPhi * iPoint._sinPhi - _cosPhi * iPoint._cosPhi * cosDif));
        }
        /// <summary>
        /// Rotates a three dimensional point about z axis.
        /// </summary>
        /// <param name="xIn">Effective X Axis for rotation.</param>
        /// <param name="yIn">Effective Y Axis for rotation.</param>
        /// <param name="zIn">Effective Z Axis for rotation.</param>
        /// <param name="cosAngle">Cosine of angle to rotate by.</param>
        /// <param name="sinAngle">Sine of angle to rotate by.</param>
        /// <returns>Position of rotated point.</returns>
        private static double[] AxisRotation(double xIn, double yIn, double zIn, double cosAngle, double sinAngle)
        {
            double[] rotatedPoint = new double[3];
            rotatedPoint[0] = xIn * cosAngle - yIn * sinAngle;
            rotatedPoint[1] = xIn * sinAngle + yIn * cosAngle;
            rotatedPoint[2] = zIn;
            return rotatedPoint;
        }
        /// <summary>
        /// Rotates a point about 3 axis in 3D Cartesian Coordinates.
        /// </summary>
        /// <param name="pointIn">Point to be rotated.</param>
        /// <param name="angle">Amount to rotate about axis.</param>
        /// <returns>Point rotated about 3 axis.</returns>
        private static double[] CartesianRotation(double[] pointIn, double[] angle)
        {
            double[] pointOut = new double[3];
            double cosAngleX = Math.Cos(angle[0]);
            double sinAngleX = Math.Sin(angle[0]);
            double cosAngleY = Math.Cos(angle[1]);
            double sinAngleY = Math.Sin(angle[1]);

            pointOut = AxisRotation(pointOut[1], pointOut[2], pointOut[0], cosAngleX, sinAngleX);
            pointOut = AxisRotation(pointOut[0], pointOut[2], pointOut[1], cosAngleY, sinAngleY);
            pointOut = AxisRotation(pointOut[0], pointOut[1], pointOut[2], Math.Cos(angle[2]), Math.Sin(angle[2]));
            pointOut = AxisRotation(pointOut[0], pointOut[2], pointOut[1], cosAngleY, -1 * sinAngleY);
            pointOut = AxisRotation(pointOut[1], pointOut[2], pointOut[0], cosAngleX, -1 * sinAngleX);

            return pointOut;
        }
        /// <summary>
        /// Rotates a point about 3 axis in 3D Spherical Coordinates.
        /// </summary>
        /// <param name="thetaIn">Input theta angle.</param>
        /// <param name="cosPhiIn">Cosine of input phi angle.</param>
        /// <param name="sinPhiIn">Sine of input phi angle.</param>
        /// <param name="angle">Angle vector to be rotated.</param>
        /// <returns>Theta and phi of rotated angle. Theta (0,2Pi], Phi (-Pi/2,Pi/2].</returns>
        private static double[] SphericalRotation(double thetaIn, double cosPhiIn, double sinPhiIn, double[] angle)
        {
            double[] pointCartesian = new double[3];

            pointCartesian[0] = Math.Cos(thetaIn) * cosPhiIn;
            pointCartesian[1] = Math.Sin(thetaIn) * cosPhiIn;
            pointCartesian[2] = sinPhiIn;

            pointCartesian = CartesianRotation(pointCartesian, angle);

            double[] pointSpherical = new double[2];

            pointSpherical[0] = (Math.Atan2(-1 * pointCartesian[1], -1 * pointCartesian[0]) + Math.PI);
            pointSpherical[1] = Math.Asin(pointCartesian[2]);
            return pointSpherical;
        }
        /// <summary>
        /// Calculates new position of point for a given rotation.
        /// </summary>
        /// <param name="angle">Three dimensional angle for rotation, given in radians.</param>
        /// <returns>New position given as X,Y.</returns>
        public int[] Transform(double[] angle)
        {
            double[] doubleOutput = SphericalRotation(_theta, _cosPhi, _sinPhi, angle);
            
            int[] outCoords = new int[2];
            outCoords[0] = (int)(Math.Round(doubleOutput[0] / _angleSize));
            if (outCoords[0] == 2*_mapSize)
            {
                outCoords[0] = 0;
            }
            outCoords[1] = (int)(Math.Round((doubleOutput[1] - _phiShift) / _angleSize));
            if (outCoords[1] == _mapSize)
            {
                outCoords[1] = 0;
            }
            return outCoords;
        }
        /// <summary>
        /// Compares two points for sorting, sorted by y first, then x.
        /// </summary>
        /// <param name="obj">Second point to compare.</param>
        /// <returns>-1 if input is after this object, 0 if the positions are equal, 
        /// and 1 if input is before this object.</returns>
        public int CompareTo(Object obj)
        {
            if (obj == null)
            {
                return 1;
            }
            if (obj is Point otherPoint)
            {
                if (this.YCoord > otherPoint.YCoord)
                {
                    return 1;
                }
                if (this.YCoord < otherPoint.YCoord)
                {
                    return -1;
                }
                if (this.XCoord > otherPoint.XCoord)
                {
                    return 1;
                }
                if (this.XCoord < otherPoint.XCoord)
                {
                    return -1;
                }
                return 0;
            }
            else
            {
                throw new ArgumentException("Object is not a Point.");
            }
        }
    }
}