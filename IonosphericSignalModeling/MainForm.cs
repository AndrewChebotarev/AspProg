namespace IonosphericSignalModeling
{
    public partial class MainForm : Form
    {
        private VideoPulse? videoPulse;
        private ImpulseResponse? impulseResponse;
        private Convolution? convolution;
        private AdaptiveFilter? adaptiveFilter;

        public MainForm()
        {
            InitializeComponent();
            InitializeClasses();
        }

        private void InitializeClasses()
        {
            videoPulse = new(this);
            impulseResponse = new(this);
            convolution = new(this, videoPulse, impulseResponse);
            adaptiveFilter = new(this, videoPulse, convolution);
        }

        private void VideoPulseSymbolTimeNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            videoPulse?.SetNewSymbolTime(videoPulseSymbolTimeNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void VideoPulseSymbolCountNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            videoPulse?.SetNewSymbolCount(videoPulseSymbolCountNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void VideoPulseAmplitudeNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            videoPulse?.SetNewAmplitude(videoPulseAmplitudeNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void ImpulseResponseMeanNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewMean(impulseResponseMeanNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void ImpulseResponseStandartDeviationNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewStandartDeviation(impulseResponseStandartDeviationNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void ImpulseResponseAmplitudeNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewAmplitude(impulseResponseAmplitudeNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void ImpulseResponseTimeNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewTime(impulseResponseTimeNumericUpDown.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void AdaptiveSilterStepNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            adaptiveFilter?.SetNewStepSize(adaptiveSilterStepNumericUpDown.Value);
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void AdaptiveSilterValueNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            adaptiveFilter?.SetNewFilterOrder(adaptiveSilterValueNumericUpDown.Value);
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox1.Checked)
            {
                checkBox2.Checked = false;
                checkBox3.Checked = false;
                groupBox5.Visible = false;
                ImpulseResponseMeanNumericUpDown_ValueChanged(new(), new());
            }
        }

        private void checkBox2_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox2.Checked)
            {
                checkBox1.Checked = false;
                checkBox3.Checked = false;
                groupBox5.Visible = true;
                impulseResponseMeanNumericUpDown2_ValueChanged(new(), new());
            }
        }

        private void checkBox3_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox3.Checked)
            {
                checkBox1.Checked = false;
                checkBox2.Checked = false;
                groupBox5.Visible = false;
            }
        }

        private void impulseResponseMeanNumericUpDown2_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewMean2(impulseResponseMeanNumericUpDown2.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void impulseResponseStandartDeviationNumericUpDown2_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewStandartDeviation2(impulseResponseStandartDeviationNumericUpDown2.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void impulseResponseAmplitudeNumericUpDown2_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewAmplitude2(impulseResponseAmplitudeNumericUpDown2.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }

        private void impulseResponseTimeNumericUpDown2_ValueChanged(object sender, EventArgs e)
        {
            impulseResponse?.SetNewTime2(impulseResponseTimeNumericUpDown2.Value);
            convolution?.SetNewConvolution();
            adaptiveFilter?.SetNewAdaptiveFilterResponse();
        }
    }
}
