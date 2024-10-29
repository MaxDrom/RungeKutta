using System.Numerics;
namespace Runge_Kutta;

public class RK4<TField, TSpace> : IODESolver<TField, TSpace>
 where TField : IFloatingPoint<TField>
 where TSpace : ILinearSpace<TSpace, TField> 
{
    private Func<TField, TSpace, TSpace> _f;
    public RK4(Func<TField, TSpace, TSpace> f)
    {
        _f = f;
    }
    public TSpace Step(TField t, TSpace x0, TField h)
    {
        var two = TField.One + TField.One;
        var k1 = _f(t, x0);
        var k2 = _f(t+h/two, x0+h/two*k1);
        var k3 = _f(t+h/two, x0+h/two*k2);
        var k4 = _f(t+h, x0 + h*k3);

        return x0 + h/(two+two+two)*(k1+two*k2+two*k3+k4);
    }
}