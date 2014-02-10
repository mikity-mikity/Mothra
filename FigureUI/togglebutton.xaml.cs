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
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace mikity.visualize
{
    /// <summary>
    /// flag.xaml の相互作用ロジック
    /// </summary>
    public partial class togglebutton : UserControl
    {
        public togglebutton()
        {
            InitializeComponent();
        }
        public string Text
        {
            set
            {
                this.button1.Content = value;
            }
            get
            {
                return (string)this.button1.Content;
            }
        }
        public bool value
        {
            get
            {
                return (bool)this.button1.IsChecked;
            }
        }
        public static implicit operator bool(mikity.visualize.togglebutton b)
        {
            return b.value;
        }
    }
}
