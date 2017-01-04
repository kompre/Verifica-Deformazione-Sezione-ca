clearvars
clc
%%
Mx = @(q,l,x) -q*x.^2/2 + q*l/2*x - q*l^2/12;
q = 100;
l = 3;
Med = Mx(q,l,0);

fi_b = 10:2:24;
sezione = rettangolo(1000,300,0,0,1,1000);
d = [40;260];
area = @(fi) 4*pi/4*fi^2;

def_not.ecu = 3.5e-3;
def_not.ec2 = 2.0e-3;
def_not.ec3 = 1.75e-3;
def_not.esu = 10e-3;
def_not.eyd = 450/1.15/210e3;

mat.cls = derivaCaratteristicheCA(25,30);
mat.steel = derivaCaratteristicheAcciaio;

fcd = mat.cls.f_cd;
fyd = mat.steel.fyd;

ris = struct('fi_b',0,'fi_a',0,'Mrd',0,'ratio',0,'lx',0,'la',0,'lt',0,'ltot',0,'peso_base',0,'peso_add',0,'peso_tot',0);
ris = repmat(ris, length(fi_b), 1);

for i = 1:length(fi_b)
    Mrd = 0;
    j = 0;
    lx = 0;
    while Mrd < abs(Med)
        fi_a = sign(j) * (fi_b(i) + 2*(j-1));
        A.tot = [area(10) ; area(fi_b(i)) + area(fi_a)];
        A.base = [area(10) ; area(fi_b(i))];
        dicotomico('x','[N, M]', 'calcoloNM(x,sezione,d,A.tot,def_not,fcd,fyd,''elastica'');', 0, 0, 300, 12)
        Mrd = M*1e-6;
        r = abs(Med)/Mrd;
        if r <= 1
            ris(i).fi_b = fi_b(i);
            ris(i).fi_a = fi_a;
            ris(i).Mrd = Mrd;
            ris(i).ratio = r;
            if fi_a > 0
                dicotomico('x','[N, M]', 'calcoloNM(x,sezione,d,A.base,def_not,fcd,fyd,''elastica'');', 0, 0, 300, 12);
                Mrd_base = M*1e-6;
                dicotomico('lx', 'Mlx', 'Mx(q,l,lx);', -Mrd_base, 0, l, 12);
                ris(i).lx = lx; % lunghezza minima necessaria 
                ris(i).la = 0.7*60*fi_a*1e-3;   % lunghezza di ancoraggio
                ris(i).lt = 0.9*max(d)/2*1e-3;   % lunghezza di traslazione
                ris(i).ltot = ris(i).lx + ris(i).la + ris(i).lt;
            end
            ris(i).peso_base = A.base(2)*1e-6 * l * 7850;
            ris(i).peso_add = 2*(area(fi_a)*1e-6 * ris(i).ltot * 7850);
            ris(i).peso_tot = ris(i).peso_base + ris(i).peso_add;
            
        end
        j = j+1;
    end
    if or(fi_a == 0, fi_b(i) == fi_a)
        break
    end
end
