namespace IonosphericSignalModeling.VideoSignalModeling._ImpulseResponse
{
    public class ImpulseResponse
    {
        private MainForm mainForm;
        private PlotModel? plotModel;
        private LineSeries? lineSeries;

        private decimal impulseResponseMean = 0.30M;
        private decimal impulseResponseStandartDeviation = 0.10M;
        private decimal impulseResponseAmplitude = 1.6M;
        private decimal impulseResponseTime = 0.7M;

        private decimal impulseResponseMean2 = 0.15M;
        private decimal impulseResponseStandardDeviation2 = 0.05M;
        private decimal impulseResponseAmplitude2 = 0.5M;
        private decimal impulseResponseTime2 = 0.5M;

        public ImpulseResponse(MainForm mainForm)
        {
            this.mainForm = mainForm;

            InitializePlot();
            InitializePoints();
        }

        public void SetNewMean(decimal newMean)
        {
            impulseResponseMean = newMean;

            if (mainForm.checkBox1.Checked)
            {
                InitializeGraph(1);
            }
            else
            {
                InitializeGraph(2);
            }
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewStandartDeviation(decimal newStandartDeviation)
        {
            impulseResponseStandartDeviation = newStandartDeviation;

            if (mainForm.checkBox1.Checked)
            {
                InitializeGraph(1);
            }
            else
            {
                InitializeGraph(2);
            }
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewAmplitude(decimal newAmplitude)
        {
            impulseResponseAmplitude = newAmplitude;

            if (mainForm.checkBox1.Checked)
            {
                InitializeGraph(1);
            }
            else
            {
                InitializeGraph(2);
            }
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewTime(decimal newTime)
        {
            impulseResponseTime = newTime;

            if (mainForm.checkBox1.Checked)
            {
                InitializeGraph(1);
            }
            else
            {
                InitializeGraph(2);
            }
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewMean2(decimal newMean2)
        {
            impulseResponseMean2 = newMean2;

            InitializeGraph(2);
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewStandartDeviation2(decimal newStandartDeviation2)
        {
            impulseResponseStandardDeviation2 = newStandartDeviation2;

            InitializeGraph(2);
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewAmplitude2(decimal newAmplitude2)
        {
            impulseResponseAmplitude2 = newAmplitude2;

            InitializeGraph(2);
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        public void SetNewTime2(decimal newTime2)
        {
            impulseResponseTime2 = newTime2;

            InitializeGraph(2);
            mainForm.impulseResponsePlotView.InvalidatePlot(true);
        }

        private void InitializePlot()
        {
            plotModel = new PlotModel();

            LinearAxis xAxis = GetXAxis();
            LinearAxis yAxis = GetYAxis();

            plotModel.Axes.Add(xAxis);
            plotModel.Axes.Add(yAxis);

            lineSeries = new LineSeries { Title = "Impulse Response" };

            plotModel.Series.Add(lineSeries);

            mainForm.impulseResponsePlotView.Model = plotModel;
        }

        private void InitializeGraph(int id)
        {
            if (id == 1)
            {
                InitializePoints();
            }
            else if (id == 2)
            {
                InitializePoints2();
            }
            else
            {

            }
        }

        private void InitializePoints()
        {
            ClearPoints();

            for (double time = 0; time <= (double)impulseResponseTime; time += 0.01)
            {
                double amplitude = GaussianFunction(time, (double)impulseResponseMean, (double)impulseResponseStandartDeviation) * (double)impulseResponseAmplitude;
                lineSeries?.Points.Add(new DataPoint(time, amplitude));
            }
        }

        private void InitializePoints2()
        {
            ClearPoints();

            // Добавляем точки для первой гауссовской кривой
            for (double time = 0; time <= (double)impulseResponseTime; time += 0.01)
            {
                double amplitude1 = GaussianFunction(time, (double)impulseResponseMean, (double)impulseResponseStandartDeviation) * (double)impulseResponseAmplitude;
                lineSeries?.Points.Add(new DataPoint(time, amplitude1));
            }

            // Добавляем точки для второй гауссовской кривой, начиная с конца первой
            double startTimeForSecondCurve = (double)impulseResponseTime + 0.01; // Начинаем с небольшого отступа
            double shiftedMean2 = startTimeForSecondCurve + (double)impulseResponseMean2;

            for (double time = startTimeForSecondCurve; time <= startTimeForSecondCurve + (double)impulseResponseTime2; time += 0.01)
            {
                double amplitude2 = GaussianFunction(time, shiftedMean2, (double)impulseResponseStandardDeviation2) * (double)impulseResponseAmplitude2;
                lineSeries?.Points.Add(new DataPoint(time, amplitude2));
            }
        }

        private double GaussianFunction(double amplitude, double mean, double standartDeviation) 
            => (1 / (standartDeviation * Math.Sqrt(2 * Math.PI))) * Math.Exp(-0.5 * Math.Pow((amplitude - mean) / standartDeviation, 2));

        private void ClearPoints() => lineSeries?.Points.Clear();

        private LinearAxis GetXAxis()
        {
            return new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Время (мс)",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.Gray,
                MinorGridlineColor = OxyColors.Gray
            };
        }

        private LinearAxis GetYAxis()
        {
            return new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Амплитуда",
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                MajorGridlineColor = OxyColors.Gray,
                MinorGridlineColor = OxyColors.Gray
            };
        }

        public List<DataPoint> GetInterpolatedPoints(int numPoints)
        {
            var interpolatedPoints = new List<DataPoint>();

            if (lineSeries == null || lineSeries.Points.Count < 2)
                return interpolatedPoints;

            for (int i = 0; i < lineSeries.Points.Count - 1; i++)
            {
                var start = lineSeries.Points[i];
                var end = lineSeries.Points[i + 1];

                interpolatedPoints.Add(start);

                for (int j = 1; j < 1; j++)
                {
                    double t = (double)j / numPoints;
                    double interpolatedX = start.X + (end.X - start.X) * t;
                    double interpolatedY = start.Y + (end.Y - start.Y) * t;
                    interpolatedPoints.Add(new DataPoint(interpolatedX, interpolatedY));
                }
            }

            interpolatedPoints.Add(lineSeries.Points.Last());
            return interpolatedPoints;
        }
    }
}