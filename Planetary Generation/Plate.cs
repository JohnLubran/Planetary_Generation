using System.Collections.Generic;

namespace Planetary_Generation
{
    /// <summary>
    /// Contains data and functions for plates.
    /// </summary>
    public class Plate
    {
        /// <summary>
        /// Magnitude of rotation per time, in radians per time unit.
        /// </summary>
        public double Speed;
        /// <summary>
        /// Two dimensional vector for rotation direction, in radians.
        /// </summary>
        public double[] Direction;
        /// <summary>
        /// Collection of points within this plate.
        /// </summary>
        public List<Point> Value;
        /// <summary>
        /// Allocates space for <see cref="Direction"/>
        /// </summary>
        public Plate()
        {
            Direction = new double[2];
            Value = new List<Point>();
        }
        /// <summary>
        /// Copies input plate points to this plate.
        /// </summary>
        /// <param name="cPlate">Input plate.</param>
        public void Copy(Plate cPlate)
        {
            Value.Clear();
            for (int i = 0; i < cPlate.Value.Count; i++)
            {
                Value.Add(new Point(cPlate.Value[i].XCoord, cPlate.Value[i].YCoord, cPlate.Value[i].Height));
            }
        }
        /// <summary>
        /// Transforms each point in plate, see <see cref="Point.Transform"/>.
        /// </summary>
        /// <param name="timeStep">Scaling factor for how much to rotate.</param>
        public void Slide(double timeStep)
        {
            int[] coords = new int[2];
            for (int i = 0; i < Value.Count; i++)
            {
                double[] angle = new double[3] { Direction[0], Direction[1], timeStep * Speed };
                coords = Value[i].Transform(angle);
                Value[i] = new Point(coords[0], coords[1], Value[i].Height);
            }
        }
        /// <summary>
        /// Sorts points and removes duplicate, averaging for every pair of duplicates with non-zero heights.
        /// </summary>
        public void CondensePoints()
        {
            Value.Sort();
            int index = 0;
            while (true)
            {
                if (Value[index].CompareTo(Value[index + 1]) == 0)
                {
                    if (Value[index].Height == 0 || Value[index + 1].Height == 0)
                    {
                        Value[index].Height += Value[index + 1].Height;
                    }
                    else
                    {
                        Value[index].Height = (Value[index].Height + Value[index + 1].Height) / 2;
                    }
                    Value.RemoveAt(index + 1);
                }
                else
                {
                    index++;
                }
                if (index >= Value.Count - 1)
                {
                    break;
                }
            }
        }
    }
}