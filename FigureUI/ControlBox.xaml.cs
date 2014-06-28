using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace Mothra.UI
{
    /// <summary>
    /// FigureUI.xaml の相互作用ロジック
    /// </summary>
    public partial class ControlBox : Window
    {
        public System.Collections.Generic.List<mikity.visualize.slider> listSlider;
        public void setCompute(Action func)
        {
            Compute.function=func;
        }
        public ControlBox()
        {
            InitializeComponent();
            listSlider = new System.Collections.Generic.List<mikity.visualize.slider>();
            this.Hide();
        }
        public void clearSliders()
        {
            this.stackPanel1.Children.Clear();
            listSlider.Clear();
        }
        public mikity.visualize.slider addSlider(int min, int step, int max,int val)
        {
            mikity.visualize.slider s=addSlider(min, step, max, val, "");
            return s;
        }
        public mikity.visualize.slider addSlider(int min, int step, int max,int val,string text)
        {
            mikity.visualize.slider s = new mikity.visualize.slider(min, step, max, val,text);
            this.stackPanel1.Children.Add(s);
            s.getSlider.ValueChanged +=new RoutedPropertyChangedEventHandler<double>(this.slider_ValueChanged);
            listSlider.Add(s);
            return s;
        }

        private void slider_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            ((mikity.visualize.slider)(from p in listSlider where (System.Object.ReferenceEquals(p.getSlider,sender))select p).First()).update();
        }

        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            e.Cancel = true;
            this.Hide();
        }
    }
}
