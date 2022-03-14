% ------------------------------------------------------------------------
% Project getting-higher/
% ------------------------------------------------------------------------
% fitnessGraphs / version 0.1 
% Matlab script to build and plot the fitness networks.
% Tested with Matlab R2018b (9.5.0.944444)
% ------------------------------------------------------------------------
% Created:  2020-08-05 
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
% IMPORTAN: The raw data comes from the simulations performed with nk_walk. 
% Before to run the present code, please verify the paths where 
% this program and your data are placed. 
% ------------------------------------------------------------------------



function inDegFitnessGraph(N,T,K,samples)
Delta_K=1;
set_of_vertices = 2^N;
while K<N
    for epi_type=1:2
        mut_op=2;
        if epi_type==1
            epi='adj';
        elseif epi_type==2
            epi='rnd';
        end
        for l=1:mut_op
            if l == 1
                mutOp = 'P';
            elseif l==2
                mutOp = 'M';
            end
            %------------------------------------------------------------------------
            G = digraph;
            %------------------------------------------------------------------------
            % Construction of the nodes genomic set
            for j = 0 : set_of_vertices-1
                s(j+1,:) = dec2bin(j,N);
                G = addnode(G,s(j+1,:));
            end
            %------------------------------------------------------------------------
            % Main loop for the analysis
            for i=1:samples
                fileID = fopen(sprintf('%s_N%dK%dT%d%s_%d', mutOp,N,K,T,epi,i),'r');
                data = textscan(fileID,'%s %f %d %d %s');
                fclose(fileID);
                f = data{2};
                genome = bin2dec(data{1});
                for k = 1 : T
                    g(k,:) = dec2bin(genome(k),N);
                end
                %--------------------------
                % assignation of fitness values per node
                for k = 1 : set_of_vertices
                    for l = 1 : T
                        if g(l,:) == s(k,:)
                            node_fitness(k) = f(l);
                        end
                    end
                end
                %--------------------------
                % construction of the fitness paths edges
                for l = 1 : T-1
                    if f(l) < f(l+1)
                        for m = 1 : set_of_vertices
                            if g(l,:) == s(m,:)
                                G = addedge(G,s(m,:),g(l+1,:));
                            end
                        end
                    end
                end
            end
            %--------------------------
            % Construction of the fitness graph
            H = simplify(G);
            h = plot(H,'LineWidth',1.5,'EdgeColor','b','NodeFontSize',12,'ArrowPosition',0.99);
            layout(h,'layered','Direction','right')
            h.EdgeAlpha = 0.35;
            wb_i = centrality(H,'outdegree');
            for ii = 1 : set_of_vertices
                if wb_i(ii) == 0
                    s(ii,:);
                    highlight(h,ii,'NodeLabel', 'r');
                end
            end
            wb_o = centrality(H,'indegree');
            for ii = 1 : set_of_vertices
                if wb_o(ii) == 0
                    s(ii,:);
                    highlight(h,ii,'NodeLabel',[ 0.9100 0.4100 0.1700]);
                end
            end
            wfit = node_fitness;
            h.NodeCData = wfit;
            colormap jet;
            c = colorbar;
            axis off
            axis square
            %--------------------------
            maxFitness = max(node_fitness);
            minFitness = min(node_fitness);
            for i = 1 : 2^N
                highlight(h,i,'MarkerSize', 20*node_fitness(i));
                if node_fitness(i) == maxFitness
                    end_node = i;
                end
                if node_fitness(i) == minFitness
                    start_node = i;
                end
            end
            %--------------------------
            p = shortestpath(H,start_node,end_node);
            highlight(h,p,'EdgeColor','r','LineWidth',2.0);            
            %title([ 'N=',num2str(N),' ','K=',num2str(K),' ',mutOp,' ',' ',epi,''], 'Interpreter','latex')
            %--------------------------
            % Verification
            q1 = ismultigraph(G);
            q2 = ismultigraph(H);
            %------------------------
            % nodes degree information
            ideg = indegree(H);
            odeg = outdegree(H);
            maxInDeg = max(ideg);
            maxOutDeg = max(odeg);
            mean_indeg = mean(ideg);
            mean_outdeg = mean(odeg);
            inCentrality = centrality(H,'indegree');
            outCentrality = centrality(H,'outdegree');
            %------------------------
            %fig = gcf;
            %saveas(gcf,'test','jpeg')
            %print('ScreenSizeFigure','-dpng','-r0')
            %title(['N=',num2str(N),' ','K=',num2str(K),' ',mutOp,' ',' ',epi,''])
            saveas(gcf,sprintf('inDegFitnessGraph_%s_N%dK%dT%d%s.png',mutOp,N,K,T,epi))
            %saveas(gcf,sprintf('inDegFitnessGraph_N%dT%d.png', N,T))
        end
    end
    K=K+Delta_K;
end
end

% counts the number of nodes with selfloops
%nodes_with_selfloops = sum(diag(Adj_mutations))

%--------------------------
% Some alternatives options:
%layout(h,'force','UseGravity',true)
%layout(h,'force','UseGravity','on')
%layout(h,'force');
%layout(h,'circle','Center',1);
%layout(h,'layered');
%layout(h,'subspace');
%layout(h,'subspace','Dimension',5)
%layout(h,'force3','Iterations',100)
%wbc = node_fitness;
%colorbar
%colorbar('southoutside');
%--------------------------
% Bye!
%--------------------------------------------------------------------------
