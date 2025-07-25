namespace IonosphericSignalModeling.VideoSignalModeling._VideoPulse
{
    public class VideoPulse
    {
        private MainForm mainForm;
        private PlotModel? plotModel;
        private LineSeries? lineSeries;

        private decimal videoPulseSymbolTime = 0.4M;
        private decimal videoPulseSymbolCount = 3;
        private decimal vidoePulseAmplitude = 1;

        public VideoPulse(MainForm mainForm)
        {
            this.mainForm = mainForm;

            InitializePlot();
            InitializePoints();
        }

        public void SetNewSymbolTime(decimal newSymbolTime)
        {
            videoPulseSymbolTime = newSymbolTime;

            mainForm.videoPulsePlotView.InvalidatePlot(true);
            InitializePoints();
        }

        public void SetNewSymbolCount(decimal newSymbolCount)
        {
            videoPulseSymbolCount = newSymbolCount;

            mainForm.videoPulsePlotView.InvalidatePlot(true);
            InitializePoints();
        }

        public void SetNewAmplitude(decimal newAmplitude)
        {
            vidoePulseAmplitude = newAmplitude;

            mainForm.videoPulsePlotView.InvalidatePlot(true);
            InitializePoints();
        }

        private void InitializePlot()
        {
            plotModel = new PlotModel();

            LinearAxis xAxis = GetXAxis();
            LinearAxis yAxis = GetYAxis();

            plotModel.Axes.Add(xAxis);
            plotModel.Axes.Add(yAxis);

            lineSeries = new LineSeries { Title = "Video Pulse" };

            plotModel.Series.Add(lineSeries);

            mainForm.videoPulsePlotView.Model = plotModel;
        }

        private void InitializePoints()
        {
            ClearPoints();

            double currentAmplitude = (double)vidoePulseAmplitude;
            double symbolDuration = (double)videoPulseSymbolTime;

            for (double time = 0; time < (double)videoPulseSymbolCount * symbolDuration; time += 0.01)
            {
                // Переключаем амплитуду между 0 и 1
                if (Math.Abs(time % symbolDuration) < 0.01)
                    currentAmplitude = (currentAmplitude == 0) ? 1 : 0;

                lineSeries?.Points.Add(new DataPoint(time, currentAmplitude));
            }

            vidoePulseAmplitude = mainForm.videoPulseAmplitudeNumericUpDown.Value;
        }

        private void ClearPoints()
        {
            lineSeries?.Points.Clear();
        }

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

                for (int j = 1; j < numPoints; j++)
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
