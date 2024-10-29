using System.Numerics;

namespace Runge_Kutta;

public class Polynom<TField>
    where TField : IFloatingPoint<TField>
{

    
    private Vector<TField> _coefs;
    public Vector<TField> Coefs => _coefs;
    public Polynom(TField[] coefs)
    {
        var ccoefs = new TField[coefs.Length];
        coefs.CopyTo(ccoefs, 0);
        _coefs = new Vector<TField>(ccoefs);
    }

    public Polynom(Vector<TField> coefs)
    {
        _coefs = coefs;
    }

    public TField Calculate(TField x)
    {
        var b = _coefs[^1];
        for (int i = _coefs.Length - 2; i>=0; i--)
            b = _coefs[i] + b*x; 

        return b;
    }

    public Polynom<TField> DivideByBinom(TField x0)
    {
        var result = new TField[_coefs.Length-1];
        var b = _coefs[^1];
        for (int i = _coefs.Length - 2; i>=0; i--)
        {
            result[i] = b;
            b = _coefs[i] + b*x0; 
        }

        return new Polynom<TField>(result);
    }

    public List<TField> GetRoots(int maxsteps = 1000)
    {
        
        var result = new List<TField>();
        if(_coefs.Length == 1)
        {
            
            return result;
        }
        if(_coefs[0] == TField.Zero)
        {
            result.Add(TField.Zero);
            result.AddRange(DivideByBinom(TField.Zero).GetRoots(maxsteps));
            return result;
        }

        if( _coefs.SubGroup(1, 2).All(x=>x==TField.Zero))
        {
            var tmp_ = new Polynom<TField>(_coefs.SubGroup(0,2).ToArray()).GetRoots(maxsteps);

            foreach(var rSquared in tmp_)
            {
                if(rSquared < TField.Zero)
                    continue;
                var r = MyMath.Root(rSquared, 2);
                result.Add(r);
                result.Add(-r);
            }

            return result;
        }

        var x = TField.Zero;
        TField dx;
        var df = DbyDx();
        var step = 0;
        do
        {
            dx = Calculate(x)/df.Calculate(x);;
            x -= dx;
            step++;
        }
        while (TField.CreateTruncating(TField.Abs(dx)) != TField.Zero && step<maxsteps);
        
        result.Add(x);
        result.AddRange(DivideByBinom(x).GetRoots(maxsteps));

        return result;
    }
    
    public Polynom<TField> DbyDx()
    {
        var result = new TField[_coefs.Length -1];
        TField ord = TField.One;
        for(var i = 0; i< result.Length; i++)
        {
            
            result[i] = ord*_coefs[i+1];
            ord++;
        }

        return new Polynom<TField>(result);
    }

    public static Polynom<TField> operator+(Polynom<TField> a, Polynom<TField> b)
    {
        var len = Math.Max(a._coefs.Length, b._coefs.Length);
        var a_tmp = new TField[len];
        var b_tmp = new TField[len];
        Array.Clear(a_tmp, 0, len);
        Array.Clear(b_tmp, 0, len);
        for(int i = 0; i< a._coefs.Length; i++)
            a_tmp[i] = a[i];

        for(int i = 0; i< b._coefs.Length; i++)
            b_tmp[i] = b[i];
        return new Polynom<TField>(new Vector<TField>(a_tmp)+ new Vector<TField>(b_tmp));
    }

    public static Polynom<TField> operator-(Polynom<TField> a, Polynom<TField> b)
    {
        var len = Math.Max(a._coefs.Length, b._coefs.Length);
        var a_tmp = new TField[len];
        var b_tmp = new TField[len];
        Array.Clear(a_tmp, 0, len);
        Array.Clear(b_tmp, 0, len);
        for(int i = 0; i< a._coefs.Length; i++)
            a_tmp[i] = a[i];

        for(int i = 0; i< b._coefs.Length; i++)
            b_tmp[i] = b[i];
        return new Polynom<TField>(new Vector<TField>(a_tmp) - new Vector<TField>(b_tmp));
    }

    public static Polynom<TField> operator*(Polynom<TField> a, Polynom<TField> b)
    {
        var result = new TField[a._coefs.Length + b._coefs.Length-1];
        var len = Math.Max(a._coefs.Length, b._coefs.Length);
        var a_tmp = new TField[len];
        var b_tmp = new TField[len];

        Array.Clear(a_tmp, 0, len);
        Array.Clear(b_tmp, 0, len);

        for(int i = 0; i< a._coefs.Length; i++)
            a_tmp[i] = a[i];

        for(int i = 0; i< b._coefs.Length; i++)
            b_tmp[i] = b[i];

        for (int i = 0; i< result.Length; i++)
        {
            result[i] = TField.Zero;
            var start = Math.Max(0, i-len +1);
            for(int j = start; (j< len) && (i-j>=0); j++)
            {
                result[i] += a_tmp[j]*b_tmp[i-j];
            }
        }
        return new Polynom<TField>(result);
    }

    public static Polynom<TField> operator*(TField a, Polynom<TField> b)
    {
        return new Polynom<TField>(new TField[]{a})*b;
    }

    public static Polynom<TField> operator*(Polynom<TField> b, TField a)
    {
        return new Polynom<TField>(new TField[]{a})*b;
    }

    public TField this[int i]
    {
        get => _coefs[i];
    }
}

internal static class Extensions
{
    public static IEnumerable<T> SubGroup<T>(this IEnumerable<T> col, int start, int step)
    {
        var en = col.Skip(start).GetEnumerator();
        var i = 0;
        while(en.MoveNext())
        {
            if(i%step == 0)
                yield return en.Current;
            i++;
        }
    }
}