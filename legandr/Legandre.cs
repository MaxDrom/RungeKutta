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
}