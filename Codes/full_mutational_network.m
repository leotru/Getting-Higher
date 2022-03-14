% ------------------------------------------------------------------------
% Project getting-higher/
% ------------------------------------------------------------------------
% full_mutational_network.m / version 0.1 
% Matlab script to build and plot the mutational networks.
% Tested with Matlab R2018b (9.5.0.944444)
% ------------------------------------------------------------------------
% Created:  2021-06-04 
% by Leonardo Trujillo (leonardo.trujillo@inria.fr)
% You may use, share, or modify this file freely
% ------------------------------------------------------------------------
% GLOBAL NOMENCLATURE:
% Compile nk_walk.cpp and run the simulation as 
% ./nk_walk  -n $N -k $K -a $A -e $epi -snk $nk_seed -swlk $wlk_seed 
%            -t $T -m $M
% where, 
% N:........Genome length.
% K:........Number of epistatic interactions (0 <= K < N).
% A:........Alphabet size (default = 2).
% epi:......Type of epistatic interactions ('ADJ' (default) or 'RND').
% nk_seed:..PRNG seed value for the landscape (default = -1).
% wlk_seed:.PRNG seed value for the random walk (default = -1).
% T:........Number of steps in the random walk (default = 10000).
% M:........Mutation type (float) (0.0 <= M <= 1.0)," 
% ..........such that M = 0.0: point mutations only (default)," 
% ....................M = 1.0: inversions only."
% ------------------------------------------------------------------------
clear;
clf;

N=8;
set_of_vertices = 2^N;

for k=1:set_of_vertices
z = (dec2bin(k-1,N));
    for j=1:N
        x(k,j) = str2double(z(j));
    end
    r = 1;
    for p=1:N
        % point mutations q=p:p / inversions q=1:N
        for q=1:N 
            y(p,q,:) = mutate(x(k,:),p,q,N);
            s(k,r,:) = strjoin(string(x(k,:)));
            t(k,r,:) = strjoin(string(y(p,q,:)));
            h(k,r) = hamming_distance(x,y(p,q,:),N); 
            r = r + 1;
        end
    end
end
mutations_graph = simplify(graph(s,t));

D = degree(mutations_graph)
dist = distances(mutations_graph);
mean_large_shortest_path = mean(max(dist));
Adj = adjacency(mutations_graph);
AdjMatrix = full(Adj);
e=eig(Adj);
c=centrality(mutations_graph,'closeness');
mg = plot(mutations_graph,'MarkerSize',10.0,'LineWidth',1.0,'NodeFontSize',12,'EdgeColor',[0.2 0.2 0.2],'NodeColor','black');
%------------------------------------------------------------------------

mg.NodeCData = D;

colormap cool;
cb = colorbar;
cb.Label.String = 'Nodes degrees';
cb.FontSize = 12;
caxis([min(D), max(D)]);

axis 'off'
axis 'square'

%--------------------------------------------------------------------------
% This function implements the Algortihm 1: Mutate(x,i,j,N).
% The procedure executes the conjugation and permutation in a single step.
% Implemented in such a way as to avoid zero-valued indices that 
% Matlab arrays do not recognize.
%--------------------------------------------------------------------------
function mutation = mutate(genome,i,j,N)
mut_genome = genome;
l = i;
while 1
    mut_temp = 1 - genome(i);
    mut_genome(i) = 1 - genome(j);
    mut_genome(j) = mut_temp;
    if j==l
        break
    end
    i = mod(i+1,N);
    if i == 0
        i=N;
    end
    j = mod(j-1,N);
    if j==0
        j=N;
    end
end
mutation = mut_genome;
end
%--------------------------------------------------------------------------
% Hamming distance
%--------------------------------------------------------------------------
function h = hamming_distance(x, y,N)
genome = x;
mut_genome = y;
h = 0;
for i=1:N
    hd = (genome(i) - mut_genome(i))^2;
    h = h + hd;
end
end

    

%--------------------------------------------------------------------------
% Bye!
%--------------------------------------------------------------------------
