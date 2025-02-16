    public static GFBiPolynomial interpolate(RSInterpolationPointSet ipSet, int dy, int k, int dy_min) throws KoetterVardyException {
        try {
            if (dy_min > dy) {
                throw new KoetterVardyException("Interpolation fails: t_min must be < t.", 4);
            }
            GaloisField field = ipSet.field();
            GFBiPolynomial oneBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{field.oneElement()}}, field);
            System.out.println("oneBiPoly: " + oneBiPoly);
            GFBiPolynomial xBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{field.zeroElement()}, {field.oneElement()}}, field);
            System.out.println("xBiPoly: " + xBiPoly);
            GFBiPolynomial[] Q = new GFBiPolynomial[dy + 1];
            int[] Wdeg = new int[dy + 1];
            GaloisField.Element[] d_rs = new GaloisField.Element[dy + 1];
            int l_min = 0;
            for (int l = 0; l <= dy; ++l) {
                Q[l] = oneBiPoly.mMul(field.oneElement(), 0, l);
                System.out.println("Q: " + Q[l]);
                Wdeg[l] = l * (k - 1);
            }
            int n = ipSet.size();
            for (int i = 0; i < n; ++i) {
                GaloisField.Element alpha_i = ipSet.pointAt(i).x();
                GaloisField.Element beta_i = ipSet.pointAt(i).y();
                int multiplicity = ipSet.pointAt(i).m();
                for (int r = 0; r < multiplicity; ++r) {
                    for (int s = 0; s < multiplicity - r; ++s) {
                        int l;
                        for (l = 0; l <= dy; ++l) {
                            d_rs[l] = Q[l].coefficient(alpha_i, beta_i, r, s);
                        }
                        l_min = KoetterVardy.argmin_l(Wdeg, d_rs);
                        if (l_min == -1) continue;
                        for (l = 0; l <= dy; ++l) {
                            if (l == l_min) continue;
                            Q[l] = Q[l].mMul(d_rs[l_min], 0, 0).sub(Q[l_min].mMul(d_rs[l], 0, 0));
                        }
                        Q[l_min] = Q[l_min].mul(xBiPoly.mSub(alpha_i, 0, 0));
                        Wdeg[l_min] = Wdeg[l_min] + 1;
                    }
                }
            }
            l_min = dy_min;
            int Wdeg_min = Wdeg[l_min];
            for (int l = 0; l <= dy; ++l) {
                if (Wdeg[l] >= Wdeg_min || Q[l].yDegree() < dy_min) continue;
                Wdeg_min = Wdeg[l];
                l_min = l;
            }
            return Q[l_min].yDegree() >= dy_min ? Q[l_min] : null;
        }
        catch (Exception e) {
            throw new KoetterVardyException("Interpolation fails.", 4, e);
        }
    }
