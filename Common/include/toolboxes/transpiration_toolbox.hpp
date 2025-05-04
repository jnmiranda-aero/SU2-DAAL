inline su2double CubicBezier(su2double B0, su2double B1,
    su2double B2, su2double B3,
    su2double eps0, su2double eps1,
    su2double eps2, su2double eps3)
{
return B0*eps0 + B1*eps1 + B2*eps2 + B3*eps3;
}

/*  Given parameter 0≤s≤1   (arc-length position along the cubic),
return εₙ.  */
inline su2double EvaluateTranspiration(const CTranspirationProfile& P,
              su2double s)
{
/* Bernstein basis */
const su2double B0 = (1-s)*(1-s)*(1-s);
const su2double B1 = 3*s*(1-s)*(1-s);
const su2double B2 = 3*s*s*(1-s);
const su2double B3 = s*s*s;

return CubicBezier(B0,B1,B2,B3, P.eps[0],P.eps[1],P.eps[2],P.eps[3]);
}
