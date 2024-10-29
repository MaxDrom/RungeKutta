using System.Numerics;

namespace Runge_Kutta;

public class BatcherRK<TField, TSpace> : IODESolver<TField, TSpace>
 where TField : IFloatingPoint<TField>
 where TSpace : ILinearSpace<TSpace, TField>
{

    private TField[,] _a;
    private TField[] _b;
    private TField[] _c;

    private int _s;
    private Polynom<TField>[] _gridPoly;

    private Func<TField, TSpace, TSpace> _f;

    private TSpace[] _lastK;

    private TField _lastH;

    public BatcherRK(Func<TField, TSpace, TSpace> f, List<TField> roots)
    {
        _s = roots.Count;
        _f = f;

        var gridPoly = new Polynom<TField>(new TField[]{-roots[0], TField.One});
        for(int i= 1; i<_s; i++)
            gridPoly *= new Polynom<TField>(new TField[]{-roots[i], TField.One});

        _gridPoly = new Polynom<TField>[_s];
        for(int i = 0; i<_s; i++)
        {
            _gridPoly[i] = gridPoly.DivideByBinom(roots[i]);
            _gridPoly[i] = _gridPoly[i] * (TField.One/_gridPoly[i].Calculate(roots[i]));
        }

        _c = [.. roots];
        for(var i = 0; i< roots.Count; i++)
            roots[i]=(roots[i] + TField.One)/(TField.One + TField.One); 

        _a = new TField[_s,_s];
        for(var i = 0; i<_s; i++)
            for (var j = 0; j<_s; j++)
                _a[i, j] = Gamma(_s, 1,  j+1, _c[i]);
        
        _b = new TField[_s];

        for(var j = 0; j<_s; j++)
            _b[j] = Gamma(_s, 1,  j+1, TField.One);
    }

    private TField Gamma(int k, int l, int j ,TField ci)
    {
        if(k == 0)
            return MyMath.Pow(ci, l)/(dynamic)MyMath.Factorial(l);
        
        if(l<= _s - j+1 && (k == j))
        {
            return Gamma(k-1, l, j, ci);
        }
        TField lf = (dynamic) l;
        return ((ci-_c[k-1])*Gamma(k-1, l, j , ci) - lf*Gamma(k-1, l+1, j , ci))/(_c[j-1] - _c[k-1]);
    }

    public TSpace Step(TField t, TSpace x0, TField h)
    {
        TSpace[] K = new TSpace[_s];
        if(_lastK != null)
        {
            for(int i = 0; i<_s; i++)
            {
                K[i] =  TSpace.AdditiveIdentity;
                for(int j = 0; j<_s; j++)
                    K[i]+=_lastK[j]*_gridPoly[j].Calculate(TField.One+_c[i]*h/_lastH);
            }
            _lastK = K;
            K = new TSpace[_s];
        }
        else
        {
            _lastK = new TSpace[_s];
            for(int i = 0; i<_s; i++)
            {
                _lastK[i] = _f(t+_c[i]*h, x0 + _c[i]*h*_f(t, x0));
            }
        }
        var tmp = new TSpace[_s];
        for(var st = 0; ; st++)
        {
            
            for(int i = 0; i<_s; i++)
            {
                tmp[i] = TSpace.AdditiveIdentity + _lastK[i];
                K[i] = TSpace.AdditiveIdentity;
                for(int j = 0; j<_s; j++)
                    K[i]+=_lastK[j]*_a[i,j];
                K[i] = x0 + h*K[i];
            }

            for(int i = 0; i<_s; i++)
                _lastK[i] = _f(t+_c[i]*h, K[i]);
            if(TField.CreateTruncating(_lastK.Zip(tmp)
                .Select(z=> z.First -z.Second)
                .Select(z=>z.Norm()).Max()) == TField.Zero)
                break;
        }

        var result = TSpace.AdditiveIdentity;
        for(int i = 0; i<_s; i++)
            result+=_lastK[i]*_b[i];
        result = x0+result*h;
        _lastH = h;
        return result;
    }
}