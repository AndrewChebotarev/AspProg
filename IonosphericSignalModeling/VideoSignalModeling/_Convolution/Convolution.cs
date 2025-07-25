namespace IonosphericSignalModeling.VideoSignalModeling._Convolution
{
    public class Convolution
    {
        private MainForm mainForm;
        private VideoPulse videoPulse;
        private ImpulseResponse impulseResponse;

        private PlotModel? plotModel;
        private LineSeries? lineSeries;
        private LineSeries? symbolSeries;

        private decimal interpolPointsCount = 1;
        private decimal discreteness = 0.01M;

        public Convolution(MainForm mainForm, VideoPulse videoPulse, ImpulseResponse impulseResponse)
        {
            this.mainForm = mainForm;
            this.videoPulse = videoPulse;
            this.impulseResponse = impulseResponse;

            InitializePlot();
            SetNewConvolution();
        }

        public void SetNewConvolution()
        {
            List<DataPoint> videoPulsePoints =  videoPulse.GetInterpolatedPoints((int)interpolPointsCount);
            List<DataPoint> impulseResponsePoints = impulseResponse.GetInterpolatedPoints((int)interpolPointsCount);

            List<DataPoint> convolutionResult = Convolve(videoPulsePoints, impulseResponsePoints);

            DisplayConvolutionResult(convolutionResult);
            DisplaySymbolTime();
        }

        private List<DataPoint> Convolve(List<DataPoint> firstSignal, List<DataPoint> secondSignal)
        {
            List<DataPoint> resultConvolve = new();
            int n = firstSignal.Count;
            int m = secondSignal.Count;

            for (int i = 0; i < n + m - 1; i++)
            {
                double sum = 0;

                for (int j = 0; j < m; j++)
                {
                    if (i - j >= 0 && i - j < n)
                        sum += firstSignal[i - j].Y * secondSignal[j].Y;
                }

                double time = i * (double)discreteness;
                resultConvolve.Add(new DataPoint(time, sum));
            }

            return resultConvolve;
        }

        private void DisplayConvolutionResult(List<DataPoint> convolutionResult)
        {
            lineSeries?.Points.Clear();
            lineSeries?.Points.AddRange(convolutionResult);
            mainForm.convolutionPlotView.InvalidatePlot(true);
        }

        private void DisplaySymbolTime()
        {
            symbolSeries?.Points.Clear();

            double minAmplitude = lineSeries?.Points.Min(p => p.Y) ?? 0;
            double maxAmplitude = lineSeries?.Points.Max(p => p.Y) ?? 0;

            symbolSeries?.Points.Add(new DataPoint(0, minAmplitude));
            symbolSeries?.Points.Add(new DataPoint(0, maxAmplitude));

            if (plotModel == null)
                return;

            var seriesToRemove = new List<LineSeries>();

            foreach (LineSeries newSymbolSeries in plotModel.Series)
            {
                if (newSymbolSeries != symbolSeries && newSymbolSeries != lineSeries)
                {
                    seriesToRemove.Add(newSymbolSeries);
                }
            }

            foreach (var series in seriesToRemove)
            {
                plotModel.Series.Remove(series);
            }

            int index = 1;
            while (index <= mainForm.videoPulseSymbolCountNumericUpDown.Value)
            {
                LineSeries newSymbolSeries = new LineSeries { Title = "Symbol visible", Color = OxyColors.Red };

                double time = (double)(index * mainForm.videoPulseSymbolTimeNumericUpDown.Value);
                newSymbolSeries.Points.Add(new DataPoint(time, minAmplitude));
                newSymbolSeries.Points.Add(new DataPoint(time, maxAmplitude));

                plotModel.Series.Add(newSymbolSeries);
                index++;
            }

            mainForm.convolutionPlotView.InvalidatePlot(true);
        }

        private void InitializePlot()
        {
            plotModel = new();

            LinearAxis xAxis = GetXAxis();
            LinearAxis yAxis = GetYAxis();

            plotModel.Axes.Add(xAxis);
            plotModel.Axes.Add(yAxis);

            lineSeries = new LineSeries { Title = "Convolution Result" };
            symbolSeries = new LineSeries { Title = "Symbol visible", Color = OxyColors.Red };

            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(symbolSeries);

            mainForm.convolutionPlotView.Model = plotModel;
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

