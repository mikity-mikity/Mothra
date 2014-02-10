using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace mikity.visualize
{
    /// <summary>
    /// Plot.xaml の相互作用ロジック
    /// </summary>
    public partial class DT : UserControl
    {
        private double dt=0.2;
        public DT()
        {
            InitializeComponent();
            Default.IsChecked=true;
        }
        public double getDt()
        {
            return dt;
        }

        private void RadioButton_Checked(object sender, RoutedEventArgs e)
        {
            dt = double.Parse((string)((RadioButton)sender).Content);
        }
    }
}
