using System.Numerics;

namespace Runge_Kutta;

public interface ILinearSpace<TSelf, TField> : IAdditionOperators<TSelf, TSelf, TSelf>,
                                               IAdditiveIdentity<TSelf, TSelf>,
                                               IMultiplyOperators<TSelf, TField, TSelf>,
                                               IUnaryNegationOperators<TSelf, TSelf>,
                                               ISubtractionOperators<TSelf,TSelf,TSelf>,
                                               IEnumerable<TField>
    where TField : INumber<TField>
    where TSelf : ILinearSpace<TSelf, TField>
{
    static abstract TSelf operator *(TField right, TSelf left);
    TField Norm();
}