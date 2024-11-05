using System.Globalization;
using System.Numerics;

namespace Runge_Kutta;

class Program
{
    static Vector<TField> Newton<TField>(TField t, Vector<TField> x)
        where TField : IFloatingPoint<TField>
    {
        var result = new Vector<TField>(4);

        result[0] = x[2];
        result[1] = x[3];

        var r = MyMath.Root(x[0]*x[0] + x[1]*x[1], 2);
        r = r*r*r;
        result[2] = -x[0]/r;
        result[3] = -x[1]/r;

        return result;

    }

    
    static void Main(string[] args)
    {
        CultureInfo.CurrentCulture = new CultureInfo("en-US", false);
        CultureInfo.CurrentCulture.NumberFormat.CurrencyDecimalDigits = 28;
        // var roots = Legandre.Polynoms<decimal>().Skip(2)
        //             .First().GetRoots();
        // for(var i = 0; i< roots.Count; i++)
        //     roots[i]=(roots[i] + decimal.One)/(decimal.One + decimal.One); 
        // roots.Sort();
        // var brk2 = new BatcherRK<decimal, Vector<decimal>>(Newton, roots);

        // roots = Legandre.Polynoms<decimal>().Skip(3)
        //             .First().GetRoots();
        // for(var i = 0; i< roots.Count; i++)
        //     roots[i]=(roots[i] + decimal.One)/(decimal.One + decimal.One); 
        // roots.Sort();
        // var brk3 = new BatcherRK<decimal, Vector<decimal>>(Newton, roots);
        // var rk = new RK4<decimal, Vector<decimal>>(Newton);
        // var integrators = new Dictionary<string, IODESolver<decimal, Vector<decimal>>>();
        // var T = 7M;
        // integrators["rk4"] = rk;
        // integrators["brk2"] = brk2;
        // integrators["brk3"] = brk3;
        // var x0 = new Vector<decimal>([0.5000000000000000000000000000M, 0.0000000000000000000000000000M,
        // 0.0000000000000000000000000000M, 1.0000000000000000000000000000M]);
        // //0.633319203086299832332011502407
        // foreach(var h in new decimal[]{0.1M, 0.01M, 0.001M, 0.0001M})
        // foreach(var (name, integrator) in integrators)
        // {
        //     using var file = new StreamWriter($"results/{h}_{name}.dat");
        //     foreach(var (t,x) in integrator.Integrate(10*T, h, x0))
        //     {
        //         file.WriteLine(t + " " + x);
        //     }
        // }
        Legandre.GetRoots<decimal>(4);
    }
}
