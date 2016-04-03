restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
with(StringTools);
with(ArrayTools);

TEST := false;
nstate := 13;
ncontr := 4;
nvar := nstate+ncontr;
nMsize := nstate*nstate;
nNsize := nstate*ncontr;
irowMxstart := nstate+1;
irowNxstart := nstate+nMsize+1;
NEQ := nstate+nMsize+nNsize;

# State und Control Vektor

x := [r[1], r[2], r[3], q[1], q[2], q[3], q[4], v[1], v[2], v[3], omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];
ydot := Vector(NEQ);
# Rotationsmatrix
R := proc (q) options operator, arrow; Matrix(3, 3, {(1, 1) = 1.-2.*q[3]^2-2.*q[4]^2, (1, 2) = -2.*q[1]*q[4]+2.*q[2]*q[3], (1, 3) = 2.*q[1]*q[3]+2.*q[2]*q[4], (2, 1) = 2.*q[1]*q[4]+2.*q[2]*q[3], (2, 2) = 1.-2.*q[2]^2-2.*q[4]^2, (2, 3) = -2.*q[1]*q[2]+2.*q[3]*q[4], (3, 1) = -2.*q[1]*q[3]+2.*q[2]*q[4], (3, 2) = 2.*q[1]*q[2]+2.*q[3]*q[4], (3, 3) = 1.-2.*q[2]^2-2.*q[3]^2}) end proc;

# Inverse der Massenmatrix
Minv := Matrix(6, 6, {(1, 1) = 1/m, (1, 2) = 0, (1, 3) = 0, (1, 4) = 0, (1, 5) = 0, (1, 6) = 0, (2, 1) = 0, (2, 2) = 1/m, (2, 3) = 0, (2, 4) = 0, (2, 5) = 0, (2, 6) = 0, (3, 1) = 0, (3, 2) = 0, (3, 3) = 1/m, (3, 4) = 0, (3, 5) = 0, (3, 6) = 0, (4, 1) = 0, (4, 2) = 0, (4, 3) = 0, (4, 4) = 1/Iges[1], (4, 5) = 0, (4, 6) = 0, (5, 1) = 0, (5, 2) = 0, (5, 3) = 0, (5, 4) = 0, (5, 5) = 1/Iges[2], (5, 6) = 0, (6, 1) = 0, (6, 2) = 0, (6, 3) = 0, (6, 4) = 0, (6, 5) = 0, (6, 6) = 1/Iges[3]});

# Kreuzproduktmatrix [q_[1:3]x]:=Matrix(3, 3, {(1, 1) = 0, (1, 2) = -q[3], (1, 3) = q[2], (2, 1) = q[3], (2, 2) = 0, (2, 3) = -q[1], (3, 1) = -q[2], (3, 2) = q[1], (3, 3) = 0});
# 
getCrossProductMatrix := proc (q) local Omega; Omega := Matrix(3, 3, shape = antisymmetric); Omega[1, 2] := -q[3]; Omega[1, 3] := q[2]; Omega[2, 3] := -q[1]; return Omega end proc;

# q(x)p := [[[p[0],-p[1:3]^(T)],[p[1:3],p[0]I[3 x3]+[p[1:3]x]]]];
quatmultiply := proc (p, q) local pTmp, qTmp, M; pTmp := `<,>`(p[1], p[2], p[3], p[4]); qTmp := `<,>`(q[1], q[2], q[3], q[4]); M := Matrix(4); M[1] := pTmp[1]; M[1, 2 .. 4] := -Transpose(pTmp[2 .. 4]); M[2 .. 4, 1] := pTmp[2 .. 4]; M[2 .. 4, 2 .. 4] := pTmp[1]*IdentityMatrix(3)+getCrossProductMatrix(pTmp[2 .. 4]); return Multiply(M, qTmp) end proc;
# Theta(q, v, omega):=[[[m*omega x v + R(q)( )^(T) * [[[0,],[0,],[g,]]]],[omega x I[ges]omega]]];
Theta := proc (q, v, omega) local res; res := Vector(6); res[1 .. 3] := m*CrossProduct(omega, v)+m*Multiply(Transpose(R(q)), `<,>`(0, 0, g)); res[4 .. 6] := CrossProduct(omega, `<,>`(Iges[1]*omega[1], Iges[2]*omega[2], Iges[3]*omega[3])); return res end proc;
# T(omega, u):= [[[0,],[0,],[-k[T]*(&sum;)u[i]^(2),],[M*u^(2)+[[[0,],[0,],[1,]]],x omega *I[M]*(u[1]-u[2]+u[3]-u[4])]]];
T := proc (omega, u) local res, M; res := Vector(6); M := Matrix(3, 4, {(1, 1) = 0, (1, 2) = d*kT, (1, 3) = 0, (1, 4) = -d*kT, (2, 1) = -d*kT, (2, 2) = 0, (2, 3) = d*kT, (2, 4) = 0, (3, 1) = -kQ, (3, 2) = kQ, (3, 3) = -kQ, (3, 4) = kQ}); res[1 .. 3] := `<,>`(0, 0, kT*(u[1]^2+u[2]^2+u[3]^2+u[4]^2)); res[4 .. 6] := Multiply(M, `<,>`(u[1]^2, u[2]^2, u[3]^2, u[4]^2))+CrossProduct(`<,>`(0, 0, 1), omega)*IM*(u[1]-u[2]+u[3]-u[4]); return res end proc;


# 
# dot(x):= [[[R(q)*v,],[0.5*q*[[[0,],[omega,]]],],[M^(-1)*(T(omega, u) - Theta(q, v, omega)),]]];
getDotBasis := proc (x) local vx, res, tmpQ, tmpV, tmpOmega, tmpU; vx := Vector(x); res := Vector(13); tmpQ := vx[4 .. 7]; tmpV := vx[8 .. 10]; tmpOmega := vx[11 .. 13]; tmpU := vx[14 .. 17]; res[1 .. 3] := Multiply(R(tmpQ), tmpV); res[4 .. 7] := .5*quatmultiply(tmpQ, `<,>`(0, tmpOmega)); res[8 .. 13] := Multiply(Minv, T(tmpOmega, tmpU)-Theta(tmpQ, tmpV, tmpOmega)); return res end proc;

getDot := proc (x) local correctTerm, res, vx, q, lambda; correctTerm := Vector(13); res := getDotBasis(x); vx := Vector(x); q := vx[4 .. 7]; lambda := 1-q[1]^2-q[2]^2-q[3]^2-q[4]^2; correctTerm[4] := lambda*q[1]; correctTerm[5] := lambda*q[2]; correctTerm[6] := lambda*q[3]; correctTerm[7] := lambda*q[4]; res := res+correctTerm; return res end proc;

# Kostenfunktion
# 
getCost := proc (x) local res, vx, s, q, v, omega, u, campos, vcampos; res := 0; vx := Vector(x); s := vx[1 .. 3]; q := vx[4 .. 7]; v := vx[8 .. 10]; omega := vx[11 .. 13]; u := vx[14 .. 17]; campos := [campos[1], campos[2], campos[3]]; vcampos := Vector(campos); res := [(1/2)*alpha*Norm(u, 2)^2+(1/2)*beta*Norm(s-vcampos, 2)^2+(1/2)*Gamma*Norm(q, 2)^2+(1/2)*kappa*Norm(vx[8 .. 13])^2]; return Vector(res) end proc;
costFunc := getCost(x);
costDFuncTmp := Transpose(Jacobian(costFunc, x));
costDFunc := Vector(nvar);
for i to nvar do costDFunc[i] := simplify(costDFuncTmp[i, 1]) end do;
costDDFunc := Hessian(costFunc[1], x);
costDDQFunc := Hessian(costFunc[1], x[1 .. 13]);
costDDRFunc := Hessian(costFunc[1], x[14 .. 17]);
# 
# 
#  Ausführen von getDot(x)

dot := getDot(x); dot := simplify(dot, 'symbolic'); dotMatrix := convert(dot, Matrix); dimDotMatrix := Dimension(dotMatrix);
# Berechne Jacobimatrix
J := Jacobian(dot, x); J := simplify(J, 'symbolic'); dimJ := Dimension(J);

# Berechne Hessematrix
H := Array(1 .. dimDotMatrix[1], 1 .. nvar, 1 .. nvar); dimH := ArrayTools[Size](H); for i to dimDotMatrix[1] do H[i] := Hessian(dot[i], x) end do;

Min := Matrix(nstate, nstate);
for i to nstate do for j to nstate do Min[i, j] := IMth(y, i, j) end do end do;

Nin := Matrix(nstate, ncontr);
for i to nstate do for j to ncontr do Nin[i, j] := INth(y, i, j) end do end do;

# Generate Funktion
Jx := J[1 .. nstate, 1 .. nstate];
Ju := J[1 .. nstate, nstate+1 .. nvar];
Mout := Multiply(Jx, Min);
Nout := Multiply(Jx, Nin)+Ju;
ydot(1 .. nstate) := dot;
ydot(nstate+1 .. nstate+nMsize) := Reshape(Mout, [nMsize, 1]);
ydot(nstate+nMsize+1 .. ()) := Reshape(Nout, [nNsize, 1]);

# Generate Jac
Jac := Matrix(NEQ, NEQ);
calc := proc (H, Msp) local k, l, m, tmp, Mspout; Mspout := Matrix(nstate, nstate); for k to nstate do for m to nstate do tmp := 0; for l to nstate do tmp := tmp+H[k, l, m]*Msp[l] end do; Mspout[k, m] := tmp end do end do; return Mspout end proc;
getNx1 := proc (H, l) local k, m, Nx; Nx := Matrix(nstate, nstate); for k to nstate do for m to nstate do Nx(k, m) := H[k, l, m] end do end do; return Nx end proc;
for i to nstate do startindexJx := (i-1)*nstate+1; endindexJx := i*nstate; startindexM := nstate+1+(i-1)*nstate; endindexM := i*nstate+nstate; startindexN := irowNxstart+(i-1)*nstate; endindexN := i*nstate+irowNxstart-1; Jac[startindexJx .. endindexJx, startindexJx .. endindexJx] := Jx; Msp := Min[1 .. nstate, i]; Mx := calc(H, Msp); Jac[startindexM .. endindexM, 1 .. nstate] := Mx; if i <= ncontr then Nsp := Nin[1 .. nstate, i]; Nx0 := calc(H, Nsp); Nx1 := getNx1(H, i+nstate); Nx := Nx0+Nx1; Jac[startindexN .. endindexN, 1 .. nstate] := Nx end if end do;
for i from nstate+1 to nvar+1 do startindexJx := (i-1)*nstate+1; endindexJx := i*nstate; Jac[startindexJx .. endindexJx, startindexJx .. endindexJx] := Jx end do;
vec := Vector(NEQ);
for i to NEQ do vec[i] := s[i] end do;

# 
# Matlabtestdaten
matlabData := "/Users/Bene/Google Drive/Documents/Studium/Master/Semester3/CaseStudies/Privat/realtime-optimization/MATLAB/tests/realtimevsfmincon";
currentdir(matlabData);
currentdir();
matlabFile := FileTools:-JoinPath(["Data.mat"], base = datadir); TestData := ImportMatrix(matlabFile, source = MATLAB, output = matrices);

NULL;
SUBNEQ := nstate;
ydot := simplify(ydot, 'symbolic', size);
Jac := simplify(Jac, 'symbolic', size);
NULL;
NULL;

TestData := TestData(1 .. NEQ);
# Temporary Ordner finden
tmpDir := TemporaryDirectory();
# Temporary Ordner als aktuellen Ordner setzen.
currentdir(tmpDir);
currentdir();
# 
configDim := Matrix([nvar, dimDotMatrix, dimJ, dimH]);
strToFileCostDD := "";

for i to nvar do for j to nvar do val := costDDFunc[i][j]; if val <> 0 then strInit := cat("cs_entry(costDD, zero_ind(", i, "), zero_ind(", j, "), "); tmp := C(val, defaulttype = float, output = string); lentmp := length(tmp); ind := SearchText("=", tmp); tmp := SubString(tmp, ind+2 .. lentmp-2); strToFileCostDD := cat(strToFileCostDD, strInit, tmp, ");\n") end if end do end do;
strToFileCostDDQ := "";
for i to nstate do for j to nstate do val := costDDQFunc[i][j]; if val <> 0 then strInit := cat("cs_entry(costDDQ, zero_ind(", i, "), zero_ind(", j, "), "); tmp := C(val, defaulttype = float, output = string); lentmp := length(tmp); ind := SearchText("=", tmp); tmp := SubString(tmp, ind+2 .. lentmp-2); strToFileCostDDQ := cat(strToFileCostDDQ, strInit, tmp, ");\n") end if end do end do;

strToFileCostDDR := "";
for i to ncontr do for j to ncontr do val := costDDRFunc[i][j]; if val <> 0 then strInit := cat("cs_entry(costDDR, zero_ind(", i, "), zero_ind(", j, "), "); tmp := C(val, defaulttype = float, output = string); lentmp := length(tmp); ind := SearchText("=", tmp); tmp := SubString(tmp, ind+2 .. lentmp-2); strToFileCostDDR := cat(strToFileCostDDR, strInit, tmp, ");\n") end if end do end do;

strToFile := "";
if TEST = false then for i to NEQ do for j to NEQ do val := Jac[i][j]; if val <> 0 then strInit := cat("IJth(pdOut, ", i, ", ", j, ")="); tmp := C(val, defaulttype = float, output = string); lentmp := length(tmp); ind := SearchText("=", tmp); tmp := SubString(tmp, ind+2 .. lentmp); strToFile := cat(strToFile, strInit, tmp) end if end do end do end if;

if TEST = false then C(ydot, defaulttype = float, resultname = ydotOut, output = tmpRTOptFunction); Text[WriteFile]("tmpRTOptJacobi", strToFile); C(costFunc[1], defaulttype = float, resultname = cost, output = tmpRTOptCost); C(costDFunc, defaulttype = float, resultname = costD, output = tmpRTOptCostD); Text[WriteFile]("tmpRTOptCostDD", strToFileCostDD); Text[WriteFile]("tmpRTOptCostDDQ", strToFileCostDDQ); Text[WriteFile]("tmpRTOptCostDDR", strToFileCostDDR) else  end if;

NULL;

