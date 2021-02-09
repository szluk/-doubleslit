function F = doubleslit(pgo1_in, pgo2_in)
% Double slit experiment based on Richard P. Feynman description
% in "QED The Strange Theory of Light and Matter".
% Illustrates the relation between distinguishability and interference. 
% INPUT:
% probabilities of opening the gaps
% pgo=1 <=> the gap is always open (always not observed); this is default value
% pgo=0 <=> the gap is always closed (always observed)
% OUTPUT:
% P - probability of detecting a particle on the detector board 
% (c) Szymon £ukaszyk
% email: szymon@patent.pl
% History
% 26.05.2020 1st published version (errors presumed)
% 09.02.2021 2nd published version

%Double slit configuration
% (x, y) - coordinates
% S      - source defining the level line (x)
% D      - detector board
% g1, g2 - top (right) and bottom (left) slit widths
% h1, h2 - height from the level line to respectively the bottom edge of the top gap and the top edge of the bottom gap 
% L1     - length from the source S to the board with gaps
% L2     - length from the board with gaps to the detector board
% dlt    - numerical resolution
% lmb    - wavelength (red EMR wavelength of 750 nm is assumed) 

% Note that we now assume zero thickness of the board with gaps. An upgrade should take it into account.
%  
%|\y          |              D,y
%|            g1,y1          | 
%|            |h1            |
%S------------|--------------|\x
%|            |h2            |
%|            g2,y2          |
%|            |              |
%|     L1     |      L2      |

switch nargin
    case 2
        pgo1 = pgo1_in;
        pgo2 = pgo2_in;
    case 1
        pgo1 = pgo1_in;
        pgo2 = 1;
    otherwise
        pgo1 = 1;
        pgo2 = 1;
end

lmb = 0.75*10^-6; %[m] (=750 nm) infrared
dlt = 0.1*10^-6;

L1 = 5000*10^-6; % = 5 mm
L2 = 5000*10^-6; % = 5 mm
h1 = 50*10^-6;   % = 50 micrometr (human hair)
h2 = 50*10^-6;   % = 50 micrometr (human hair)
g1 = 10*10^-6;   % = 10 micrometr
g2 = 10*10^-6;   % = 10 micrometr

H1 =  h1; HG1 = h1+g1;    % top gap coordinates
H2 = -h2; HG2 = -(h2+g2); % bottom gap coordinates
y1 = H1 :dlt:HG1;         % variable coordinate of the top gap
y2 = HG2:dlt:H2;          % variable coordinate of the bottom gap
y = -5*h1:dlt:5*h2;       % variable coordinate of the detector

% calculate paths between the source and the detector board for each detector position
for k=1:length(y)
  d11 = sqrt(L1^2 + y1.^2);         % from the source to the top gap
  d21 = sqrt(L2^2 + (y(k)-y1).^2);  % from the top gap to the detector board
  d1 = d11+d21;

  d12 = sqrt(L1^2 + y2.^2);         % from the source to the bottom gap
  d22 = sqrt(L2^2 + (y(k)-y2).^2);  % from the bottom gap to the detector board
  d2 = d12+d22;

  % calculate phasors
  d1 = 2*pi*d1/lmb;
  d2 = 2*pi*d2/lmb;

  % calculate arrows s coordinates sx, sy; one arrow corresponds to one path
  sx1 = cos(d1);      
  sy1 = sin(d1);

  sx2 = cos(d2);
  sy2 = sin(d2);

  dx=0;
  dy=0;
  ns=0; % number of arrows

  b1 = pgo1; a1 = b1-1; % range of uniform probability distribution of the top gap detector
  b2 = pgo2; a2 = b2-1; % range of uniform probability distribution of the bottom gap detector

  for l=1:length(sx1)
    go1 = a1 + (b1-a1).*rand(1);
    if go1 > 0
      dx = dx + sx1(l);
      dy = dy + sy1(l);  
      ns = ns+1;
    end
  end
  for l=1:length(sx2)
    go2 = a2 + (b2-a2).*rand(1);
    if go2 > 0 
      dx = dx + sx2(l);
      dy = dy + sy2(l);    
      ns = ns+1;    
    end
  end

  len = sqrt( ( dx./ns ).^2 + (dy./ns).^2);
  P(k)= len.^2;  % probability of detecting a particle on the detector board
end

figure
F=plot(y,P)
x_lab=sprintf('Screen width [m]\n(L_1=L_2=5 mm; h_1=h_2=50 µm; g_1=g_2=10 µm; \\lambda=750 nm [infrared])');

xlabel(x_lab)
ylabel('Detection probability');
tit_lab = sprintf('Gap observation probabilities: p(left)=%3.2f, p(right)=%3.2f', 1-pgo2_in, 1-pgo1_in);
title(tit_lab);

axis([min(y) max(y) 0 1]);

line([H1  H1 ], [0 1],'LineStyle',':')
line([HG1 HG1], [0 1],'LineStyle',':')

line([H2  H2 ], [0 1],'LineStyle',':')
line([HG2 HG2], [0 1],'LineStyle',':')

%line([0 1*10^-4], [-1 1],'LineStyle',':')

return
