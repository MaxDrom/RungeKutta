using System.Numerics;

namespace Runge_Kutta;

public interface IODESolver<TField, TSpace>
 where TField : IFloatingPoint<TField>
 where TSpace : ILinearSpace<TSpace, TField>
{
    TSpace Step(TField t, TSpace x0, TField h);
}

public static class IODESolverExtensions
{
    public static IEnumerable<(TField, TSpace)> Integrate<TField, TSpace>(this IODESolver<TField, TSpace> integrator, 
            TField T, TField h, TSpace x0)
        where TField : IFloatingPoint<TField>
        where TSpace : ILinearSpace<TSpace, TField>
    {
        TField t = TField.Zero;
        var x = x0;
        do
        {
            yield return (t, x);
            x = integrator.Step(t, x, h);
            t+=h;
        }
        while (t <= T);
    }
}