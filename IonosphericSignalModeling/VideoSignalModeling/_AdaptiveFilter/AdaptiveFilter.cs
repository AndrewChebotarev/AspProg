namespace IonosphericSignalModeling.VideoSignalModeling._AdaptiveFilter
{
    public class AdaptiveFilter
    {
        private MainForm mainForm;
        private VideoPulse videoPulse;
        private Convolution convolution;

        private PlotModel? plotModel;
        private LineSeries? lineSeries;

        private decimal interpolPointsCount = 100;
        private decimal discreteness = 0.01M;

        private double stepSize = 0.5;
        private int filterOrder = 1;
        private double[] coefficients;

        public AdaptiveFilter(MainForm mainForm, VideoPulse videoPulse, Convolution convolution)
        {
            this.mainForm = mainForm;
            this.videoPulse = videoPulse;
            this.convolution = convolution;

            coefficients = new double[filterOrder];
            InitializePlot();
            SetNewAdaptiveFilterResponse();
        }

        public void SetNewStepSize(decimal newStepSize)
        {
            this.stepSize = (double)newStepSize;
            SetNewAdaptiveFilterResponse();
        }

        public void SetNewFilterOrder(decimal newFilterOrder)
        {
            try
            {
                this.filterOrder = (int)newFilterOrder;
                coefficients = new double[filterOrder];
                SetNewAdaptiveFilterResponse();
            }
            catch
            {

            }
        }

        public void SetNewAdaptiveFilterResponse()
        {
            List<DataPoint> videoPulsePoints = videoPulse.GetInterpolatedPoints((int)interpolPointsCount);
            List<DataPoint> convolutionResultPoints = convolution.GetInterpolatedPoints((int)interpolPointsCount);

            List<DataPoint> adaptiveFilterResponse = ApplyAdaptiveFilter(convolutionResultPoints, videoPulsePoints);

            DisplayImpulseResponse(adaptiveFilterResponse);
        }

        private List<DataPoint> ApplyAdaptiveFilter(List<DataPoint> convolutionResultPoints, List<DataPoint> videoPulsePoints)
        {
            List<DataPoint> resultImpulseResponse = new();
            int n = convolutionResultPoints.Count;

            // Паддинг видео импульса
            List<DataPoint> paddedVideoPulse = new List<DataPoint>(videoPulsePoints);
            for (int i = videoPulsePoints.Count; i < n; i++)
            {
                paddedVideoPulse.Add(new DataPoint(i * (double)discreteness, 0));
            }

            // Инициализация коэффициентов
            coefficients = new double[filterOrder];
            for (int i = 0; i < filterOrder; i++)
            {
                coefficients[i] = 0;
            }

            for (int i = 0; i < n; i++)
            {
                double output = 0;

                // Вычисление выхода фильтра
                for (int j = 0; j < filterOrder; j++)
                {
                    if (i - j >= 0)
                    {
                        output += coefficients[j] * convolutionResultPoints[i - j].Y;
                    }
                }

                // Ошибка между желаемым и текущим выходом
                double desired = paddedVideoPulse[i].Y;
                double error = desired - output;

                // Обновление коэффициентов без нормализации
                for (int j = 0; j < filterOrder; j++)
                {
                    if (i - j >= 0)
                    {
                        coefficients[j] += stepSize * error * convolutionResultPoints[i - j].Y;
                    }
                }

                // Сохранение результата
                double time = i * (double)discreteness;
                resultImpulseResponse.Add(new DataPoint(time, output));
            }

            return resultImpulseResponse;
        }

        private void DisplayImpulseResponse(List<DataPoint> impulseResponseResult)
        {
            lineSeries?.Points.Clear();
            lineSeries?.Points.AddRange(impulseResponseResult);
            mainForm.adaptiveFilterPlotView.InvalidatePlot(true);
        }

        private void InitializePlot()
        {
            plotModel = new();

            LinearAxis xAxis = GetXAxis();
            LinearAxis yAxis = GetYAxis();

            plotModel.Axes.Add(xAxis);
            plotModel.Axes.Add(yAxis);

            lineSeries = new LineSeries { Title = "Impulse Response from Adaptive Filter" };

            plotModel.Series.Add(lineSeries);

            mainForm.adaptiveFilterPlotView.Model = plotModel;
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
    }
}