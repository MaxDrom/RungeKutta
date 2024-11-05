using System.Numerics;
using Runge_Kutta;

public static class Legandre
{
    public static IEnumerable<Polynom<TField>> Polynoms<TField>()
        where TField : IFloatingPoint<TField>
    {
        var p1 = new Polynom<TField>(new TField[]{TField.One});
        yield return p1;
        var p2 = new Polynom<TField>(new TField[]{TField.Zero, TField.One});
        yield return p2;
        var order = TField.One;
        var x = new Polynom<TField>(new TField[]{TField.Zero, TField.One});
        var two = TField.One + TField.One;
        while( true)
        {
            var p3 = (two*order+TField.One)/(TField.One + order)*x*p2 - order/(TField.One + order)*p1;
            order++;
            p1 = p2;
            p2 = p3;
            yield return p3;
        }
    }

    public static IEnumerable<TField> Polynoms<TField>(TField x)
        where TField : IFloatingPoint<TField>
    {
        var p1 = TField.One;
        yield return p1;
        var p2 = x;
        yield return p2;
        var order = TField.One;
        var two = TField.One + TField.One;
        while( true)
        {
            var p3 = (two*order+TField.One)/(TField.One + order)*x*p2 - order/(TField.One + order)*p1;
            order++;
            p1 = p2;
            p2 = p3;
            yield return p3;
        }
    } 
    public static List<TField> GetRoots<TField>(int order)
        where TField : IFloatingPoint<TField>
    {
        TField forder = (dynamic) (order/2);
        var results = new List<TField>();
        var two = TField.One + TField.One;
        var four = two+two;
        if(order%2 == 1)
            results.Add(TField.Zero);
        for(TField i = TField.One; i<=forder; i++)
        {
            var d = TField.Pi*(four*i-TField.One)/(four*forder*two + two);
            var x =MyMath.Sin(d);
            
            x = MyMath.Root(TField.One-x*x,2);

            for(int st = 0; st<100; st++)
                x -= Polynoms(x).Skip(order).First()/DPolynoms(x).Skip(order).First();

            results.Add(x);
            results.Add(-x);
        }

        foreach(var r in results)
            Console.WriteLine(r+" " + Polynoms(r).Skip(order).First());
        return results;
        
    }
    public static IEnumerable<TField> DPolynoms<TField>(TField x)
        where TField : IFloatingPoint<TField>
    {
        yield return TField.Zero;
        var order = TField.One;
        var denom = TField.One - x*x;
        foreach(var (pp, p) in Polynoms(x).Zip(Polynoms(x).Skip(1)))
        {
            yield return order/denom*(pp - x*p);
            order++;
        }
    }
}