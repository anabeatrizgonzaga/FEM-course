import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import sys

def main():
    print("=== Programa MEF - Metodo dos Elementos Finitos ===")
    fname_in = input("Planilha de entrada (.xlsx): ").strip()
    fname_out = input("Arquivo de saida (VTK): ").strip()

    print("\nLendo dados da planilha...")
    nodes, elements, materials, bcs, forces, ndf, ndm, nen = read_excel_data(fname_in)
    nnode = len(nodes)
    neq = nnode * ndf

    print("Montando a matriz de rigidez global...")
    K_global = lil_matrix((neq, neq))
    F_global = np.zeros(neq)

    for node_idx, force_vec in forces.items():
        for j in range(ndf):
            F_global[(node_idx-1)*ndf + j] += force_vec[j]

    for el in elements:
        ix = el['ix']
        mat = materials[el['mat']]
        xl = np.array([nodes[n] for n in ix])
        Ke = get_element_stiffness(xl, mat, ndf)
        for i in range(len(ix)):
            for j in range(len(ix)):
                for a in range(ndf):
                    for b in range(ndf):
                        row = (ix[i]-1)*ndf + a
                        col = (ix[j]-1)*ndf + b
                        K_global[row, col] += Ke[i*ndf+a, j*ndf+b]

    print("Aplicando condicoes de contorno...")
    PENALTY = 1e15
    for node_idx, bc_vec in bcs.items():
        for j in range(ndf):
            if bc_vec[j] == 1:
                eq_idx = (node_idx-1)*ndf + j
                K_global[eq_idx, eq_idx] = PENALTY
                F_global[eq_idx] = 0.0

    # --- Solver direto (robusto com penalizacao alta) ---
    print("Resolvendo o sistema (solver direto)...")
    K_csc = K_global.tocsc()
    U_global = spsolve(K_csc, F_global)
    print("Sistema resolvido.")

    # --- Recuperacao de tensoes ---
    print("Calculando tensoes nodais...")
    stress = compute_nodal_stresses(nodes, elements, materials, U_global, ndf)
    vm = von_mises(stress)

    print("Gerando arquivo VTK...")
    write_vtk(fname_out, nodes, elements, U_global, ndf, stress, vm)
    print("Calculo concluido com sucesso!")

# ======================================================================
# Rotinas de Elementos Finitos (rigidez)
# ======================================================================

def get_element_stiffness(xl, mat, ndf):
    iel = mat['iel']
    if iel in [1, 2]:
        return elmt_t3(xl, mat, iel)
    elif iel in [3, 4]:
        return elmt_q4(xl, mat, iel)
    elif iel == 7:
        return elmt_q8(xl, mat, iel)
    else:
        raise ValueError(f"Elemento tipo {iel} ainda nao implementado.")

def elmt_t3(xl, mat, iel):
    E=mat['E']; nu=mat['nu']; thic=mat.get('thic',1.0)
    if iel==1:
        c=E*(1-nu)/((1+nu)*(1-2*nu))
        D=np.array([[c,c*nu/(1-nu),0],[c*nu/(1-nu),c,0],[0,0,E/(2*(1+nu))]])
    else:
        c=E/(1-nu**2); D=np.array([[c,c*nu,0],[c*nu,c,0],[0,0,E/(2*(1+nu))]])
    x1,y1=xl[0]; x2,y2=xl[1]; x3,y3=xl[2]
    area=0.5*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
    b1,c1=y2-y3,x3-x2; b2,c2=y3-y1,x1-x3; b3,c3=y1-y2,x2-x1
    B=(1.0/(2.0*area))*np.array([[b1,0,b2,0,b3,0],[0,c1,0,c2,0,c3],[c1,b1,c2,b2,c3,b3]])
    return B.T@D@B*area*thic

def elmt_q4(xl, mat, iel):
    E=mat['E']; nu=mat['nu']; thic=mat.get('thic',1.0)
    if iel==3:
        c=E*(1-nu)/((1+nu)*(1-2*nu))
        D=np.array([[c,c*nu/(1-nu),0],[c*nu/(1-nu),c,0],[0,0,E/(2*(1+nu))]])
    else:
        c=E/(1-nu**2); D=np.array([[c,c*nu,0],[c*nu,c,0],[0,0,E/(2*(1+nu))]])
    pts=[-0.577350269189626,0.577350269189626]; Ke=np.zeros((8,8))
    for r in pts:
        for s in pts:
            dhdr=np.array([(1+s),-(1+s),-(1-s),(1-s)])/4.0
            dhds=np.array([(1+r),(1-r),-(1-r),-(1+r)])/4.0
            J=np.zeros((2,2))
            J[0,0]=np.dot(dhdr,xl[:,0]); J[0,1]=np.dot(dhdr,xl[:,1])
            J[1,0]=np.dot(dhds,xl[:,0]); J[1,1]=np.dot(dhds,xl[:,1])
            detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            B=np.zeros((3,8))
            for i in range(4):
                d=invJ@np.array([dhdr[i],dhds[i]])
                B[0,2*i]=d[0]; B[1,2*i+1]=d[1]; B[2,2*i]=d[1]; B[2,2*i+1]=d[0]
            Ke+=B.T@D@B*detJ*thic
    return Ke

def elmt_q8(xl, mat, iel):
    E=mat['E']; nu=mat['nu']; thic=mat.get('thic',1.0)
    c=E/(1-nu**2); D=np.array([[c,c*nu,0],[c*nu,c,0],[0,0,E/(2*(1+nu))]])
    pts=[-np.sqrt(0.6),0.0,np.sqrt(0.6)]; wts=[5.0/9.0,8.0/9.0,5.0/9.0]
    Ke=np.zeros((16,16))
    for i in range(3):
        for j in range(3):
            r=pts[i]; s=pts[j]; weight=wts[i]*wts[j]
            dhdr=np.array([0.25*(1-s)*(2*r+s),0.25*(1-s)*(2*r-s),0.25*(1+s)*(2*r+s),
                0.25*(1+s)*(2*r-s),-r*(1-s),0.5*(1-s**2),-r*(1+s),-0.5*(1-s**2)])
            dhds=np.array([0.25*(1-r)*(2*s+r),0.25*(1+r)*(2*s-r),0.25*(1+r)*(2*s+r),
                0.25*(1-r)*(2*s-r),-0.5*(1-r**2),-s*(1+r),0.5*(1-r**2),-s*(1-r)])
            J=np.zeros((2,2))
            J[0,0]=np.dot(dhdr,xl[:,0]); J[0,1]=np.dot(dhdr,xl[:,1])
            J[1,0]=np.dot(dhds,xl[:,0]); J[1,1]=np.dot(dhds,xl[:,1])
            detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            B=np.zeros((3,16))
            for k in range(8):
                d=invJ@np.array([dhdr[k],dhds[k]])
                B[0,2*k]=d[0]; B[1,2*k+1]=d[1]; B[2,2*k]=d[1]; B[2,2*k+1]=d[0]
            Ke+=B.T@D@B*detJ*thic*weight
    return Ke

# ======================================================================
# Recuperacao de tensoes (T3, Q4, Q8)
# ======================================================================

# Coordenadas naturais dos nos, NA MESMA ORDENACAO usada em cada elemento
NODE_NATURAL = {
    3: np.array([[0,0],[0,0],[0,0]], dtype=float),               # T3: B constante
    4: np.array([[1,1],[-1,1],[-1,-1],[1,-1]], dtype=float),      # Q4 (ordem do codigo)
    8: np.array([[-1,-1],[1,-1],[1,1],[-1,1],
                 [0,-1],[1,0],[0,1],[-1,0]], dtype=float),        # Q8 (ordem do codigo)
}

def get_D(E, nu, iel):
    if iel in (1, 3):  # EPD
        c=E*(1-nu)/((1+nu)*(1-2*nu))
        return np.array([[c,c*nu/(1-nu),0],[c*nu/(1-nu),c,0],[0,0,E/(2*(1+nu))]])
    c=E/(1-nu**2)      # EPT
    return np.array([[c,c*nu,0],[c*nu,c,0],[0,0,E/(2*(1+nu))]])

def _shape_derivs(nen, r, s):
    if nen==4:
        dhdr=np.array([(1+s),-(1+s),-(1-s),(1-s)])/4.0
        dhds=np.array([(1+r),(1-r),-(1-r),-(1+r)])/4.0
    else:  # nen==8
        dhdr=np.array([0.25*(1-s)*(2*r+s),0.25*(1-s)*(2*r-s),0.25*(1+s)*(2*r+s),
            0.25*(1+s)*(2*r-s),-r*(1-s),0.5*(1-s**2),-r*(1+s),-0.5*(1-s**2)])
        dhds=np.array([0.25*(1-r)*(2*s+r),0.25*(1+r)*(2*s-r),0.25*(1+r)*(2*s+r),
            0.25*(1-r)*(2*s-r),-0.5*(1-r**2),-s*(1+r),0.5*(1-r**2),-s*(1-r)])
    return dhdr, dhds

def _bmatrix(nen, r, s, xl):
    if nen==3:
        x1,y1=xl[0]; x2,y2=xl[1]; x3,y3=xl[2]
        area=0.5*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
        b1,c1=y2-y3,x3-x2; b2,c2=y3-y1,x1-x3; b3,c3=y1-y2,x2-x1
        return (1/(2*area))*np.array([[b1,0,b2,0,b3,0],[0,c1,0,c2,0,c3],[c1,b1,c2,b2,c3,b3]])
    dhdr,dhds=_shape_derivs(nen,r,s)
    J=np.zeros((2,2))
    J[0,0]=dhdr@xl[:,0]; J[0,1]=dhdr@xl[:,1]
    J[1,0]=dhds@xl[:,0]; J[1,1]=dhds@xl[:,1]
    invJ=np.linalg.inv(J)
    B=np.zeros((3,nen*2))
    for k in range(nen):
        d=invJ@np.array([dhdr[k],dhds[k]])
        B[0,2*k]=d[0]; B[1,2*k+1]=d[1]; B[2,2*k]=d[1]; B[2,2*k+1]=d[0]
    return B

def compute_nodal_stresses(nodes, elements, materials, U, ndf):
    """sxx, syy, txy nos nos por avaliacao direta + media entre elementos vizinhos."""
    nnode=len(nodes)
    ssum=np.zeros((nnode+1,3)); cnt=np.zeros(nnode+1)
    for el in elements:
        ix=el['ix']; mat=materials[el['mat']]
        D=get_D(mat['E'],mat['nu'],mat['iel'])
        xl=np.array([nodes[n] for n in ix]); nen=len(ix)
        ue=np.zeros(nen*ndf)
        for a in range(nen):
            ue[a*ndf]=U[(ix[a]-1)*ndf]; ue[a*ndf+1]=U[(ix[a]-1)*ndf+1]
        nat=NODE_NATURAL[nen]
        for a in range(nen):
            r,s=nat[a]
            sig=D@(_bmatrix(nen,r,s,xl)@ue)
            ssum[ix[a]]+=sig; cnt[ix[a]]+=1
    stress=np.zeros((nnode+1,3))
    for n in range(1,nnode+1):
        if cnt[n]>0: stress[n]=ssum[n]/cnt[n]
    return stress

def von_mises(stress):
    sxx,syy,txy=stress[:,0],stress[:,1],stress[:,2]
    return np.sqrt(sxx**2 - sxx*syy + syy**2 + 3*txy**2)

# ======================================================================
# Leitura e Escrita
# ======================================================================

def read_excel_data(filename):
    try:
        df_nos=pd.read_excel(filename,sheet_name='Nos')
        df_elem=pd.read_excel(filename,sheet_name='Elementos')
        df_mat=pd.read_excel(filename,sheet_name='Materiais')
        df_cc=pd.read_excel(filename,sheet_name='Contorno')
        df_forcas=pd.read_excel(filename,sheet_name='Forcas')
    except Exception as e:
        sys.exit(f"Erro ao abrir a planilha: {e}")
    ndm=2; ndf=2
    if 'N8' in df_elem.columns:
        nen=8; cols=['N1','N2','N3','N4','N5','N6','N7','N8']
        print(f"-> Elemento Q8 detectado ({nen} nos por elemento).")
    else:
        nen=4; cols=['N1','N2','N3','N4']
        print(f"-> Elemento Q4/T3 detectado ({nen} nos por elemento).")
    nodes={int(r['ID']):[float(r['X']),float(r['Y'])] for _,r in df_nos.iterrows()}
    elements=[]
    for _,r in df_elem.iterrows():
        ix=[int(r[c]) for c in cols if not pd.isna(r.get(c,np.nan))]
        elements.append({'id':int(r['ID']),'ix':ix,'mat':int(r['Mat'])})
    materials={}
    for _,r in df_mat.iterrows():
        materials[int(r['ID'])]={'iel':int(r['Tipo']),'E':float(r['E']),'nu':float(r['nu']),
            'thic':float(r['Espessura']) if 'Espessura' in r and not pd.isna(r['Espessura']) else 1.0}
    bcs={int(r['ID']):[int(r['Restr_X']),int(r['Restr_Y'])] for _,r in df_cc.iterrows()}
    forces={int(r['ID']):[float(r['Fx']),float(r['Fy'])] for _,r in df_forcas.iterrows()}
    return nodes,elements,materials,bcs,forces,ndf,ndm,nen

def write_vtk(filename, nodes, elements, U, ndf, stress=None, vm=None):
    nnode=len(nodes); numel=len(elements)
    with open(filename,'w') as f:
        f.write("# vtk DataFile Version 2.0\nMEF - Python 2D\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {nnode} float\n")
        for i in range(1,nnode+1):
            f.write(f"{nodes[i][0]} {nodes[i][1]} 0.0\n")
        nen=len(elements[0]['ix'])
        f.write(f"CELLS {numel} {numel*(nen+1)}\n")
        for el in elements:
            f.write(f"{nen} " + " ".join(str(n-1) for n in el['ix']) + "\n")
        f.write(f"CELL_TYPES {numel}\n")
        cell_type={3:5,4:9,8:23}[nen]
        for _ in range(numel): f.write(f"{cell_type}\n")
        # --- Dados nodais ---
        f.write(f"POINT_DATA {nnode}\n")
        f.write("VECTORS Deslocamentos float\n")
        for i in range(1,nnode+1):
            ux=U[(i-1)*ndf]; uy=U[(i-1)*ndf+1] if ndf>1 else 0.0
            f.write(f"{ux:.8e} {uy:.8e} 0.0\n")
        if stress is not None:
            for col,nome in [(0,"Sxx"),(1,"Syy"),(2,"Txy")]:
                f.write(f"SCALARS {nome} float 1\nLOOKUP_TABLE default\n")
                for i in range(1,nnode+1):
                    f.write(f"{stress[i,col]:.8e}\n")
        if vm is not None:
            f.write("SCALARS VonMises float 1\nLOOKUP_TABLE default\n")
            for i in range(1,nnode+1):
                f.write(f"{vm[i]:.8e}\n")

if __name__ == "__main__":
    main()
