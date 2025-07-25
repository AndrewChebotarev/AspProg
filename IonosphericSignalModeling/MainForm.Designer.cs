namespace IonosphericSignalModeling
{
    partial class MainForm
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(MainForm));
            tabControl1 = new TabControl();
            tabPage1 = new TabPage();
            groupBox6 = new GroupBox();
            label6 = new Label();
            impulseResponseMeanNumericUpDown = new NumericUpDown();
            label5 = new Label();
            impulseResponseTimeNumericUpDown = new NumericUpDown();
            impulseResponseStandartDeviationNumericUpDown = new NumericUpDown();
            label7 = new Label();
            label4 = new Label();
            impulseResponseAmplitudeNumericUpDown = new NumericUpDown();
            groupBox5 = new GroupBox();
            impulseResponseTimeNumericUpDown2 = new NumericUpDown();
            label18 = new Label();
            impulseResponseAmplitudeNumericUpDown2 = new NumericUpDown();
            label19 = new Label();
            impulseResponseStandartDeviationNumericUpDown2 = new NumericUpDown();
            label20 = new Label();
            impulseResponseMeanNumericUpDown2 = new NumericUpDown();
            label21 = new Label();
            groupBox4 = new GroupBox();
            groupBox3 = new GroupBox();
            adaptiveSilterValueNumericUpDown = new NumericUpDown();
            label10 = new Label();
            adaptiveSilterStepNumericUpDown = new NumericUpDown();
            label11 = new Label();
            groupBox2 = new GroupBox();
            checkBox3 = new CheckBox();
            checkBox2 = new CheckBox();
            checkBox1 = new CheckBox();
            groupBox1 = new GroupBox();
            videoPulseAmplitudeNumericUpDown = new NumericUpDown();
            label3 = new Label();
            videoPulseSymbolCountNumericUpDown = new NumericUpDown();
            label2 = new Label();
            videoPulseSymbolTimeNumericUpDown = new NumericUpDown();
            label1 = new Label();
            adaptiveFilterPlotView = new OxyPlot.WindowsForms.PlotView();
            convolutionPlotView = new OxyPlot.WindowsForms.PlotView();
            impulseResponsePlotView = new OxyPlot.WindowsForms.PlotView();
            videoPulsePlotView = new OxyPlot.WindowsForms.PlotView();
            tabPage2 = new TabPage();
            tabControl1.SuspendLayout();
            tabPage1.SuspendLayout();
            groupBox6.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)impulseResponseMeanNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseTimeNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseStandartDeviationNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseAmplitudeNumericUpDown).BeginInit();
            groupBox5.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)impulseResponseTimeNumericUpDown2).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseAmplitudeNumericUpDown2).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseStandartDeviationNumericUpDown2).BeginInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseMeanNumericUpDown2).BeginInit();
            groupBox3.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)adaptiveSilterValueNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)adaptiveSilterStepNumericUpDown).BeginInit();
            groupBox2.SuspendLayout();
            groupBox1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)videoPulseAmplitudeNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)videoPulseSymbolCountNumericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)videoPulseSymbolTimeNumericUpDown).BeginInit();
            SuspendLayout();
            // 
            // tabControl1
            // 
            tabControl1.Controls.Add(tabPage1);
            tabControl1.Controls.Add(tabPage2);
            tabControl1.Dock = DockStyle.Fill;
            tabControl1.Location = new Point(0, 0);
            tabControl1.Name = "tabControl1";
            tabControl1.SelectedIndex = 0;
            tabControl1.Size = new Size(1924, 1055);
            tabControl1.TabIndex = 0;
            // 
            // tabPage1
            // 
            tabPage1.Controls.Add(groupBox6);
            tabPage1.Controls.Add(groupBox5);
            tabPage1.Controls.Add(groupBox4);
            tabPage1.Controls.Add(groupBox3);
            tabPage1.Controls.Add(groupBox2);
            tabPage1.Controls.Add(groupBox1);
            tabPage1.Controls.Add(adaptiveFilterPlotView);
            tabPage1.Controls.Add(convolutionPlotView);
            tabPage1.Controls.Add(impulseResponsePlotView);
            tabPage1.Controls.Add(videoPulsePlotView);
            tabPage1.Location = new Point(4, 29);
            tabPage1.Name = "tabPage1";
            tabPage1.Padding = new Padding(3);
            tabPage1.Size = new Size(1916, 1022);
            tabPage1.TabIndex = 0;
            tabPage1.Text = "Видеоимпульс";
            tabPage1.UseVisualStyleBackColor = true;
            // 
            // groupBox6
            // 
            groupBox6.Controls.Add(label6);
            groupBox6.Controls.Add(impulseResponseMeanNumericUpDown);
            groupBox6.Controls.Add(label5);
            groupBox6.Controls.Add(impulseResponseTimeNumericUpDown);
            groupBox6.Controls.Add(impulseResponseStandartDeviationNumericUpDown);
            groupBox6.Controls.Add(label7);
            groupBox6.Controls.Add(label4);
            groupBox6.Controls.Add(impulseResponseAmplitudeNumericUpDown);
            groupBox6.Location = new Point(1147, 239);
            groupBox6.Name = "groupBox6";
            groupBox6.Size = new Size(384, 221);
            groupBox6.TabIndex = 9;
            groupBox6.TabStop = false;
            groupBox6.Text = "Параметры гаус. кривой";
            // 
            // label6
            // 
            label6.AutoSize = true;
            label6.Location = new Point(21, 34);
            label6.Name = "label6";
            label6.Size = new Size(142, 20);
            label6.TabIndex = 6;
            label6.Text = "Средние значение:";
            // 
            // impulseResponseMeanNumericUpDown
            // 
            impulseResponseMeanNumericUpDown.DecimalPlaces = 2;
            impulseResponseMeanNumericUpDown.Increment = new decimal(new int[] { 5, 0, 0, 131072 });
            impulseResponseMeanNumericUpDown.Location = new Point(213, 32);
            impulseResponseMeanNumericUpDown.Minimum = new decimal(new int[] { 100, 0, 0, int.MinValue });
            impulseResponseMeanNumericUpDown.Name = "impulseResponseMeanNumericUpDown";
            impulseResponseMeanNumericUpDown.Size = new Size(150, 27);
            impulseResponseMeanNumericUpDown.TabIndex = 7;
            impulseResponseMeanNumericUpDown.Value = new decimal(new int[] { 30, 0, 0, 131072 });
            impulseResponseMeanNumericUpDown.ValueChanged += ImpulseResponseMeanNumericUpDown_ValueChanged;
            // 
            // label5
            // 
            label5.AutoSize = true;
            label5.Location = new Point(21, 67);
            label5.Name = "label5";
            label5.Size = new Size(186, 20);
            label5.TabIndex = 8;
            label5.Text = "Средквадрат отклонение:";
            // 
            // impulseResponseTimeNumericUpDown
            // 
            impulseResponseTimeNumericUpDown.DecimalPlaces = 1;
            impulseResponseTimeNumericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            impulseResponseTimeNumericUpDown.Location = new Point(213, 133);
            impulseResponseTimeNumericUpDown.Maximum = new decimal(new int[] { 10000, 0, 0, 0 });
            impulseResponseTimeNumericUpDown.Name = "impulseResponseTimeNumericUpDown";
            impulseResponseTimeNumericUpDown.Size = new Size(150, 27);
            impulseResponseTimeNumericUpDown.TabIndex = 13;
            impulseResponseTimeNumericUpDown.Value = new decimal(new int[] { 7, 0, 0, 65536 });
            impulseResponseTimeNumericUpDown.ValueChanged += ImpulseResponseTimeNumericUpDown_ValueChanged;
            // 
            // impulseResponseStandartDeviationNumericUpDown
            // 
            impulseResponseStandartDeviationNumericUpDown.DecimalPlaces = 2;
            impulseResponseStandartDeviationNumericUpDown.Increment = new decimal(new int[] { 5, 0, 0, 131072 });
            impulseResponseStandartDeviationNumericUpDown.Location = new Point(213, 67);
            impulseResponseStandartDeviationNumericUpDown.Name = "impulseResponseStandartDeviationNumericUpDown";
            impulseResponseStandartDeviationNumericUpDown.Size = new Size(150, 27);
            impulseResponseStandartDeviationNumericUpDown.TabIndex = 9;
            impulseResponseStandartDeviationNumericUpDown.Value = new decimal(new int[] { 10, 0, 0, 131072 });
            impulseResponseStandartDeviationNumericUpDown.ValueChanged += ImpulseResponseStandartDeviationNumericUpDown_ValueChanged;
            // 
            // label7
            // 
            label7.AutoSize = true;
            label7.Location = new Point(21, 135);
            label7.Name = "label7";
            label7.Size = new Size(57, 20);
            label7.TabIndex = 12;
            label7.Text = "Время:";
            // 
            // label4
            // 
            label4.AutoSize = true;
            label4.Location = new Point(21, 97);
            label4.Name = "label4";
            label4.Size = new Size(79, 20);
            label4.TabIndex = 10;
            label4.Text = "Амлитуда:";
            // 
            // impulseResponseAmplitudeNumericUpDown
            // 
            impulseResponseAmplitudeNumericUpDown.DecimalPlaces = 1;
            impulseResponseAmplitudeNumericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            impulseResponseAmplitudeNumericUpDown.Location = new Point(213, 100);
            impulseResponseAmplitudeNumericUpDown.Maximum = new decimal(new int[] { 10000, 0, 0, 0 });
            impulseResponseAmplitudeNumericUpDown.Name = "impulseResponseAmplitudeNumericUpDown";
            impulseResponseAmplitudeNumericUpDown.Size = new Size(150, 27);
            impulseResponseAmplitudeNumericUpDown.TabIndex = 11;
            impulseResponseAmplitudeNumericUpDown.Value = new decimal(new int[] { 16, 0, 0, 65536 });
            impulseResponseAmplitudeNumericUpDown.ValueChanged += ImpulseResponseAmplitudeNumericUpDown_ValueChanged;
            // 
            // groupBox5
            // 
            groupBox5.Controls.Add(impulseResponseTimeNumericUpDown2);
            groupBox5.Controls.Add(label18);
            groupBox5.Controls.Add(impulseResponseAmplitudeNumericUpDown2);
            groupBox5.Controls.Add(label19);
            groupBox5.Controls.Add(impulseResponseStandartDeviationNumericUpDown2);
            groupBox5.Controls.Add(label20);
            groupBox5.Controls.Add(impulseResponseMeanNumericUpDown2);
            groupBox5.Controls.Add(label21);
            groupBox5.Location = new Point(1537, 239);
            groupBox5.Name = "groupBox5";
            groupBox5.Size = new Size(391, 221);
            groupBox5.TabIndex = 8;
            groupBox5.TabStop = false;
            groupBox5.Text = "Параметры гаус. кривой - 2";
            groupBox5.Visible = false;
            // 
            // impulseResponseTimeNumericUpDown2
            // 
            impulseResponseTimeNumericUpDown2.DecimalPlaces = 1;
            impulseResponseTimeNumericUpDown2.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            impulseResponseTimeNumericUpDown2.Location = new Point(211, 131);
            impulseResponseTimeNumericUpDown2.Maximum = new decimal(new int[] { 10000, 0, 0, 0 });
            impulseResponseTimeNumericUpDown2.Name = "impulseResponseTimeNumericUpDown2";
            impulseResponseTimeNumericUpDown2.Size = new Size(150, 27);
            impulseResponseTimeNumericUpDown2.TabIndex = 13;
            impulseResponseTimeNumericUpDown2.Value = new decimal(new int[] { 5, 0, 0, 65536 });
            impulseResponseTimeNumericUpDown2.ValueChanged += impulseResponseTimeNumericUpDown2_ValueChanged;
            // 
            // label18
            // 
            label18.AutoSize = true;
            label18.Location = new Point(19, 133);
            label18.Name = "label18";
            label18.Size = new Size(57, 20);
            label18.TabIndex = 12;
            label18.Text = "Время:";
            // 
            // impulseResponseAmplitudeNumericUpDown2
            // 
            impulseResponseAmplitudeNumericUpDown2.DecimalPlaces = 1;
            impulseResponseAmplitudeNumericUpDown2.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            impulseResponseAmplitudeNumericUpDown2.Location = new Point(211, 98);
            impulseResponseAmplitudeNumericUpDown2.Maximum = new decimal(new int[] { 10000, 0, 0, 0 });
            impulseResponseAmplitudeNumericUpDown2.Name = "impulseResponseAmplitudeNumericUpDown2";
            impulseResponseAmplitudeNumericUpDown2.Size = new Size(150, 27);
            impulseResponseAmplitudeNumericUpDown2.TabIndex = 11;
            impulseResponseAmplitudeNumericUpDown2.Value = new decimal(new int[] { 5, 0, 0, 65536 });
            impulseResponseAmplitudeNumericUpDown2.ValueChanged += impulseResponseAmplitudeNumericUpDown2_ValueChanged;
            // 
            // label19
            // 
            label19.AutoSize = true;
            label19.Location = new Point(19, 95);
            label19.Name = "label19";
            label19.Size = new Size(79, 20);
            label19.TabIndex = 10;
            label19.Text = "Амлитуда:";
            // 
            // impulseResponseStandartDeviationNumericUpDown2
            // 
            impulseResponseStandartDeviationNumericUpDown2.DecimalPlaces = 2;
            impulseResponseStandartDeviationNumericUpDown2.Increment = new decimal(new int[] { 5, 0, 0, 131072 });
            impulseResponseStandartDeviationNumericUpDown2.Location = new Point(211, 65);
            impulseResponseStandartDeviationNumericUpDown2.Name = "impulseResponseStandartDeviationNumericUpDown2";
            impulseResponseStandartDeviationNumericUpDown2.Size = new Size(150, 27);
            impulseResponseStandartDeviationNumericUpDown2.TabIndex = 9;
            impulseResponseStandartDeviationNumericUpDown2.Value = new decimal(new int[] { 5, 0, 0, 131072 });
            impulseResponseStandartDeviationNumericUpDown2.ValueChanged += impulseResponseStandartDeviationNumericUpDown2_ValueChanged;
            // 
            // label20
            // 
            label20.AutoSize = true;
            label20.Location = new Point(19, 65);
            label20.Name = "label20";
            label20.Size = new Size(186, 20);
            label20.TabIndex = 8;
            label20.Text = "Средквадрат отклонение:";
            // 
            // impulseResponseMeanNumericUpDown2
            // 
            impulseResponseMeanNumericUpDown2.DecimalPlaces = 2;
            impulseResponseMeanNumericUpDown2.Increment = new decimal(new int[] { 5, 0, 0, 131072 });
            impulseResponseMeanNumericUpDown2.Location = new Point(211, 30);
            impulseResponseMeanNumericUpDown2.Minimum = new decimal(new int[] { 100, 0, 0, int.MinValue });
            impulseResponseMeanNumericUpDown2.Name = "impulseResponseMeanNumericUpDown2";
            impulseResponseMeanNumericUpDown2.Size = new Size(150, 27);
            impulseResponseMeanNumericUpDown2.TabIndex = 7;
            impulseResponseMeanNumericUpDown2.Value = new decimal(new int[] { 15, 0, 0, 131072 });
            impulseResponseMeanNumericUpDown2.ValueChanged += impulseResponseMeanNumericUpDown2_ValueChanged;
            // 
            // label21
            // 
            label21.AutoSize = true;
            label21.Location = new Point(19, 32);
            label21.Name = "label21";
            label21.Size = new Size(142, 20);
            label21.TabIndex = 6;
            label21.Text = "Средние значение:";
            // 
            // groupBox4
            // 
            groupBox4.Location = new Point(831, 466);
            groupBox4.Name = "groupBox4";
            groupBox4.Size = new Size(700, 221);
            groupBox4.TabIndex = 7;
            groupBox4.TabStop = false;
            groupBox4.Text = "Параметры свертки";
            // 
            // groupBox3
            // 
            groupBox3.Controls.Add(adaptiveSilterValueNumericUpDown);
            groupBox3.Controls.Add(label10);
            groupBox3.Controls.Add(adaptiveSilterStepNumericUpDown);
            groupBox3.Controls.Add(label11);
            groupBox3.Location = new Point(831, 693);
            groupBox3.Name = "groupBox3";
            groupBox3.Size = new Size(700, 221);
            groupBox3.TabIndex = 6;
            groupBox3.TabStop = false;
            groupBox3.Text = "Параметры адапитивного фильтра";
            // 
            // adaptiveSilterValueNumericUpDown
            // 
            adaptiveSilterValueNumericUpDown.Location = new Point(160, 73);
            adaptiveSilterValueNumericUpDown.Maximum = new decimal(new int[] { 1000, 0, 0, 0 });
            adaptiveSilterValueNumericUpDown.Minimum = new decimal(new int[] { 100, 0, 0, int.MinValue });
            adaptiveSilterValueNumericUpDown.Name = "adaptiveSilterValueNumericUpDown";
            adaptiveSilterValueNumericUpDown.Size = new Size(150, 27);
            adaptiveSilterValueNumericUpDown.TabIndex = 13;
            adaptiveSilterValueNumericUpDown.Value = new decimal(new int[] { 1, 0, 0, 0 });
            adaptiveSilterValueNumericUpDown.ValueChanged += AdaptiveSilterValueNumericUpDown_ValueChanged;
            // 
            // label10
            // 
            label10.AutoSize = true;
            label10.Location = new Point(21, 75);
            label10.Name = "label10";
            label10.Size = new Size(127, 20);
            label10.TabIndex = 12;
            label10.Text = "Порядок филтра:";
            // 
            // adaptiveSilterStepNumericUpDown
            // 
            adaptiveSilterStepNumericUpDown.DecimalPlaces = 7;
            adaptiveSilterStepNumericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 458752 });
            adaptiveSilterStepNumericUpDown.Location = new Point(160, 32);
            adaptiveSilterStepNumericUpDown.Minimum = new decimal(new int[] { 100, 0, 0, int.MinValue });
            adaptiveSilterStepNumericUpDown.Name = "adaptiveSilterStepNumericUpDown";
            adaptiveSilterStepNumericUpDown.Size = new Size(150, 27);
            adaptiveSilterStepNumericUpDown.TabIndex = 11;
            adaptiveSilterStepNumericUpDown.Value = new decimal(new int[] { 5, 0, 0, 65536 });
            adaptiveSilterStepNumericUpDown.ValueChanged += AdaptiveSilterStepNumericUpDown_ValueChanged;
            // 
            // label11
            // 
            label11.AutoSize = true;
            label11.Location = new Point(21, 34);
            label11.Name = "label11";
            label11.Size = new Size(111, 20);
            label11.TabIndex = 10;
            label11.Text = "Шаг обучения:";
            // 
            // groupBox2
            // 
            groupBox2.Controls.Add(checkBox3);
            groupBox2.Controls.Add(checkBox2);
            groupBox2.Controls.Add(checkBox1);
            groupBox2.Location = new Point(831, 239);
            groupBox2.Name = "groupBox2";
            groupBox2.Size = new Size(310, 221);
            groupBox2.TabIndex = 5;
            groupBox2.TabStop = false;
            groupBox2.Text = "Параметры импульсной характеристики";
            // 
            // checkBox3
            // 
            checkBox3.AutoSize = true;
            checkBox3.Location = new Point(21, 86);
            checkBox3.Name = "checkBox3";
            checkBox3.Size = new Size(176, 24);
            checkBox3.TabIndex = 16;
            checkBox3.Text = "Ионосферный канал";
            checkBox3.UseVisualStyleBackColor = true;
            checkBox3.CheckedChanged += checkBox3_CheckedChanged;
            // 
            // checkBox2
            // 
            checkBox2.AutoSize = true;
            checkBox2.Location = new Point(21, 56);
            checkBox2.Name = "checkBox2";
            checkBox2.Size = new Size(145, 24);
            checkBox2.TabIndex = 15;
            checkBox2.Text = "Гаус. кривые = 2";
            checkBox2.UseVisualStyleBackColor = true;
            checkBox2.CheckedChanged += checkBox2_CheckedChanged;
            // 
            // checkBox1
            // 
            checkBox1.AutoSize = true;
            checkBox1.Checked = true;
            checkBox1.CheckState = CheckState.Checked;
            checkBox1.Location = new Point(21, 26);
            checkBox1.Name = "checkBox1";
            checkBox1.Size = new Size(142, 24);
            checkBox1.TabIndex = 14;
            checkBox1.Text = "Гаус. кривая = 1";
            checkBox1.UseVisualStyleBackColor = true;
            checkBox1.CheckedChanged += checkBox1_CheckedChanged;
            // 
            // groupBox1
            // 
            groupBox1.Controls.Add(videoPulseAmplitudeNumericUpDown);
            groupBox1.Controls.Add(label3);
            groupBox1.Controls.Add(videoPulseSymbolCountNumericUpDown);
            groupBox1.Controls.Add(label2);
            groupBox1.Controls.Add(videoPulseSymbolTimeNumericUpDown);
            groupBox1.Controls.Add(label1);
            groupBox1.Location = new Point(831, 31);
            groupBox1.Name = "groupBox1";
            groupBox1.Size = new Size(700, 202);
            groupBox1.TabIndex = 4;
            groupBox1.TabStop = false;
            groupBox1.Text = "Параметры видеоимпульса";
            // 
            // videoPulseAmplitudeNumericUpDown
            // 
            videoPulseAmplitudeNumericUpDown.Location = new Point(231, 110);
            videoPulseAmplitudeNumericUpDown.Maximum = new decimal(new int[] { 10000, 0, 0, 0 });
            videoPulseAmplitudeNumericUpDown.Name = "videoPulseAmplitudeNumericUpDown";
            videoPulseAmplitudeNumericUpDown.Size = new Size(150, 27);
            videoPulseAmplitudeNumericUpDown.TabIndex = 5;
            videoPulseAmplitudeNumericUpDown.Value = new decimal(new int[] { 1, 0, 0, 0 });
            videoPulseAmplitudeNumericUpDown.ValueChanged += VideoPulseAmplitudeNumericUpDown_ValueChanged;
            // 
            // label3
            // 
            label3.AutoSize = true;
            label3.Location = new Point(21, 117);
            label3.Name = "label3";
            label3.Size = new Size(79, 20);
            label3.TabIndex = 4;
            label3.Text = "Амлитуда:";
            // 
            // videoPulseSymbolCountNumericUpDown
            // 
            videoPulseSymbolCountNumericUpDown.Location = new Point(231, 67);
            videoPulseSymbolCountNumericUpDown.Name = "videoPulseSymbolCountNumericUpDown";
            videoPulseSymbolCountNumericUpDown.Size = new Size(150, 27);
            videoPulseSymbolCountNumericUpDown.TabIndex = 3;
            videoPulseSymbolCountNumericUpDown.Value = new decimal(new int[] { 3, 0, 0, 0 });
            videoPulseSymbolCountNumericUpDown.ValueChanged += VideoPulseSymbolCountNumericUpDown_ValueChanged;
            // 
            // label2
            // 
            label2.AutoSize = true;
            label2.Location = new Point(21, 74);
            label2.Name = "label2";
            label2.Size = new Size(166, 20);
            label2.TabIndex = 2;
            label2.Text = "Количество символов:";
            // 
            // videoPulseSymbolTimeNumericUpDown
            // 
            videoPulseSymbolTimeNumericUpDown.DecimalPlaces = 1;
            videoPulseSymbolTimeNumericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            videoPulseSymbolTimeNumericUpDown.Location = new Point(231, 26);
            videoPulseSymbolTimeNumericUpDown.Name = "videoPulseSymbolTimeNumericUpDown";
            videoPulseSymbolTimeNumericUpDown.Size = new Size(150, 27);
            videoPulseSymbolTimeNumericUpDown.TabIndex = 1;
            videoPulseSymbolTimeNumericUpDown.Value = new decimal(new int[] { 4, 0, 0, 65536 });
            videoPulseSymbolTimeNumericUpDown.ValueChanged += VideoPulseSymbolTimeNumericUpDown_ValueChanged;
            // 
            // label1
            // 
            label1.AutoSize = true;
            label1.Location = new Point(21, 33);
            label1.Name = "label1";
            label1.Size = new Size(204, 20);
            label1.TabIndex = 0;
            label1.Text = "Длительность символа (мс):";
            // 
            // adaptiveFilterPlotView
            // 
            adaptiveFilterPlotView.Location = new Point(33, 693);
            adaptiveFilterPlotView.Name = "adaptiveFilterPlotView";
            adaptiveFilterPlotView.PanCursor = Cursors.Hand;
            adaptiveFilterPlotView.Size = new Size(792, 221);
            adaptiveFilterPlotView.TabIndex = 3;
            adaptiveFilterPlotView.Text = "plotView4";
            adaptiveFilterPlotView.ZoomHorizontalCursor = Cursors.SizeWE;
            adaptiveFilterPlotView.ZoomRectangleCursor = Cursors.SizeNWSE;
            adaptiveFilterPlotView.ZoomVerticalCursor = Cursors.SizeNS;
            // 
            // convolutionPlotView
            // 
            convolutionPlotView.Location = new Point(33, 466);
            convolutionPlotView.Name = "convolutionPlotView";
            convolutionPlotView.PanCursor = Cursors.Hand;
            convolutionPlotView.Size = new Size(792, 221);
            convolutionPlotView.TabIndex = 2;
            convolutionPlotView.Text = "plotView3";
            convolutionPlotView.ZoomHorizontalCursor = Cursors.SizeWE;
            convolutionPlotView.ZoomRectangleCursor = Cursors.SizeNWSE;
            convolutionPlotView.ZoomVerticalCursor = Cursors.SizeNS;
            // 
            // impulseResponsePlotView
            // 
            impulseResponsePlotView.Location = new Point(33, 239);
            impulseResponsePlotView.Name = "impulseResponsePlotView";
            impulseResponsePlotView.PanCursor = Cursors.Hand;
            impulseResponsePlotView.Size = new Size(792, 221);
            impulseResponsePlotView.TabIndex = 1;
            impulseResponsePlotView.Text = "plotView2";
            impulseResponsePlotView.ZoomHorizontalCursor = Cursors.SizeWE;
            impulseResponsePlotView.ZoomRectangleCursor = Cursors.SizeNWSE;
            impulseResponsePlotView.ZoomVerticalCursor = Cursors.SizeNS;
            // 
            // videoPulsePlotView
            // 
            videoPulsePlotView.Location = new Point(33, 31);
            videoPulsePlotView.Name = "videoPulsePlotView";
            videoPulsePlotView.PanCursor = Cursors.Hand;
            videoPulsePlotView.Size = new Size(792, 202);
            videoPulsePlotView.TabIndex = 0;
            videoPulsePlotView.Text = "plotView1";
            videoPulsePlotView.ZoomHorizontalCursor = Cursors.SizeWE;
            videoPulsePlotView.ZoomRectangleCursor = Cursors.SizeNWSE;
            videoPulsePlotView.ZoomVerticalCursor = Cursors.SizeNS;
            // 
            // tabPage2
            // 
            tabPage2.Location = new Point(4, 29);
            tabPage2.Name = "tabPage2";
            tabPage2.Padding = new Padding(3);
            tabPage2.Size = new Size(1916, 1022);
            tabPage2.TabIndex = 1;
            tabPage2.Text = "Другое";
            tabPage2.UseVisualStyleBackColor = true;
            // 
            // MainForm
            // 
            AutoScaleDimensions = new SizeF(8F, 20F);
            AutoScaleMode = AutoScaleMode.Font;
            ClientSize = new Size(1924, 1055);
            Controls.Add(tabControl1);
            Icon = (Icon)resources.GetObject("$this.Icon");
            Name = "MainForm";
            Text = "ISM";
            tabControl1.ResumeLayout(false);
            tabPage1.ResumeLayout(false);
            groupBox6.ResumeLayout(false);
            groupBox6.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)impulseResponseMeanNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseTimeNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseStandartDeviationNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseAmplitudeNumericUpDown).EndInit();
            groupBox5.ResumeLayout(false);
            groupBox5.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)impulseResponseTimeNumericUpDown2).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseAmplitudeNumericUpDown2).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseStandartDeviationNumericUpDown2).EndInit();
            ((System.ComponentModel.ISupportInitialize)impulseResponseMeanNumericUpDown2).EndInit();
            groupBox3.ResumeLayout(false);
            groupBox3.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)adaptiveSilterValueNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)adaptiveSilterStepNumericUpDown).EndInit();
            groupBox2.ResumeLayout(false);
            groupBox2.PerformLayout();
            groupBox1.ResumeLayout(false);
            groupBox1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)videoPulseAmplitudeNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)videoPulseSymbolCountNumericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)videoPulseSymbolTimeNumericUpDown).EndInit();
            ResumeLayout(false);
        }

        #endregion

        private TabControl tabControl1;
        private TabPage tabPage1;
        private TabPage tabPage2;
        public OxyPlot.WindowsForms.PlotView convolutionPlotView;
        public OxyPlot.WindowsForms.PlotView impulseResponsePlotView;
        public OxyPlot.WindowsForms.PlotView viplotView1;
        private GroupBox groupBox1;
        private GroupBox groupBox3;
        private GroupBox groupBox2;
        public OxyPlot.WindowsForms.PlotView videoPulsePlotView;
        public NumericUpDown videoPulseSymbolCountNumericUpDown;
        private Label label2;
        public NumericUpDown videoPulseSymbolTimeNumericUpDown;
        private Label label1;
        public NumericUpDown videoPulseAmplitudeNumericUpDown;
        private Label label3;
        public NumericUpDown impulseResponseAmplitudeNumericUpDown;
        private Label label4;
        private NumericUpDown impulseResponseStandartDeviationNumericUpDown;
        private Label label5;
        private NumericUpDown impulseResponseMeanNumericUpDown;
        private Label label6;
        public NumericUpDown impulseResponseTimeNumericUpDown;
        private Label label7;
        private GroupBox groupBox4;
        private NumericUpDown convolutionDiscNumericUpDown;
        private NumericUpDown convolutionInterpolNumericUpDown;
        public OxyPlot.WindowsForms.PlotView adaptiveFilterPlotView;
        private NumericUpDown adaptiveSilterValueNumericUpDown;
        private Label label10;
        private NumericUpDown adaptiveSilterStepNumericUpDown;
        private Label label11;
        public CheckBox checkBox3;
        public CheckBox checkBox2;
        public CheckBox checkBox1;
        private GroupBox groupBox6;
        private GroupBox groupBox5;
        public NumericUpDown impulseResponseTimeNumericUpDown2;
        private Label label18;
        public NumericUpDown impulseResponseAmplitudeNumericUpDown2;
        private Label label19;
        private NumericUpDown impulseResponseStandartDeviationNumericUpDown2;
        private Label label20;
        private NumericUpDown impulseResponseMeanNumericUpDown2;
        private Label label21;
    }
}
