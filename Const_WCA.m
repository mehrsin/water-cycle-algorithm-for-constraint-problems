function [Xmin,Fmin,Sum_Const,NFEs,Elapsed_Time]=Const_WCA(objective_function,constraints,LB,UB,nvars,Npop,Nsr,dmax,max_it)
% Water Cycle Algorithm (WCA) for solving constrained minimization problems
%% Inputs

% objective_function                    Objective function which you wish to optimize (Max,Min)
% constraints                           Imposed constraints of a given problem
% LB                                    Lower bound of a given problem
% UB                                    Upper bound of a given problem
% nvars                                 Number of design variables
% Npop                                  Population size
% Nsr                                   Number of rivers and sea
% dmax                                  Evaporation condition constant
% max_it                                Maximum number of iteration

% %  Outputs

% Xmin                                  Obtained optimum solution
% Fmin                                  Cost (fitness) of obtained optimum solution
% Sum_Const                             Summation of constraint violations
% NFEs                                  Number of function evaluations
% Elapsed_Time                          Elapsed time

%% Defualt User parameters
format long g
if (nargin<6 || isempty(Npop)), Npop=50; end
if (nargin<7 || isempty(Nsr)), Nsr=4; end
if (nargin<8 || isempty(dmax)), dmax=1e-5; end    % For constrained problems
if (nargin<9 || isempty(max_it)), max_it=1000; end

%% Create initial population and form sea,rivers, and streams
tic;
N_stream=Npop-Nsr;

epss=eps;

ind.position=[];
ind.cost=[];
ind.const=[];

pop=repmat(ind,Npop,1);

for i=1:Npop
    pop(i).position=LB+(UB-LB).*rand(1,nvars);
    pop(i).cost=objective_function(pop(i).position);
    c=constraints(pop(i).position);
    pop(i).const=sum(c(c>epss)); 
end

X_MINUS=[];
aa=[pop.const];
COST_MINUS=[pop(aa<=epss).cost];
if ~isempty(COST_MINUS)
    X_MINUS=pop(aa<=epss);
    [~,INDEX_M]=sort(COST_MINUS);
    X_MINUS=X_MINUS(INDEX_M);
end

X_PLUS=[];
SUM_C_PLUS=aa(aa>epss);
if ~isempty(SUM_C_PLUS)
    X_PLUS=pop(aa>epss);
    [~,INDEX_P]=sort(SUM_C_PLUS);
    X_PLUS=X_PLUS(INDEX_P);
end

pop=[X_MINUS;X_PLUS];
% ----------------------------- Forming the Sea ---------------------------
sea=pop(1);
% ----------------------------- Forming the Rivers ------------------------
river=pop(2:Nsr);
% ----------------------------- Forming the Streams -----------------------
stream=pop(Nsr+1:end);
% -------------- Designate streams to the rivers and sea ------------------
cs=[sea.cost;[river.cost]';stream(1).cost];

CN=cs-max(cs);

NS=round(abs(CN/sum(CN))*N_stream);
NS(end)=[];
NS=sort(NS,'descend');
% --------------------- Modification on NS --------------------------------
i=Nsr;
while sum(NS)>N_stream
    if NS(i)>1
        NS(i)=NS(i)-1;
    else
        i=i-1;
    end
end

i=1;
while sum(NS)<N_stream
    NS(i)=NS(i)+1;
end

if find(NS==0)
    index=find(NS==0);
    for i=1:size(index,1)
        while NS(index(i))==0
            NS(index(i))=NS(index(i))+round(NS(i)/6);
            NS(i)=NS(i)-round(NS(i)/6);
        end
    end
end

NS=sort(NS,'descend');
NB=NS(2:end);
%% Main Loop WCA
FF=zeros(max_it,1);
for i=1:max_it
    % Moving Streams to the Sea
    for j=1:NS(1)
        stream(j).position=stream(j).position+2.*rand.*(sea.position-stream(j).position);
        
        stream(j).position=min(stream(j).position,UB);
        stream(j).position=max(stream(j).position,LB);
        
        stream(j).cost=objective_function(stream(j).position);
        c=constraints(stream(j).position);
        stream(j).const=sum(c(c>epss));
        
        if (sea.const<=epss && stream(j).const<=epss && stream(j).cost<sea.cost) || (sea.const>epss && stream(j).const<=epss)
            new_sea=stream(j);
            stream(j)=sea;
            sea=new_sea;
        elseif sea.const>epss && stream(j).const>epss && stream(j).const<sea.const
            new_sea=stream(j);
            stream(j)=sea;
            sea=new_sea;
        end
    end
    % Moving Streams to Rivers
    for k=1:Nsr-1
        for j=1:NB(k)
            stream(j+sum(NS(1:k))).position=stream(j+sum(NS(1:k))).position+rand(1,nvars).*2.*(river(k).position-stream(j+sum(NS(1:k))).position);
            
            stream(j+sum(NS(1:k))).position=min(stream(j+sum(NS(1:k))).position,UB);
            stream(j+sum(NS(1:k))).position=max(stream(j+sum(NS(1:k))).position,LB);
            
            stream(j+sum(NS(1:k))).cost=objective_function(stream(j+sum(NS(1:k))).position);
            c=constraints(stream(j+sum(NS(1:k))).position);
            stream(j+sum(NS(1:k))).const=sum(c(c>epss));
            
            Yes=0;
            if (river(k).const<=epss && stream(j+sum(NS(1:k))).const<=epss && stream(j+sum(NS(1:k))).cost<river(k).cost) || (river(k).const>epss && stream(j+sum(NS(1:k))).const<=eps)
                new_river=stream(j+sum(NS(1:k)));
                stream(j+sum(NS(1:k)))=river(k);
                river(k)=new_river;
                Yes=1;
            elseif river(k).const>epss && stream(j+sum(NS(1:k))).const>epss && stream(j+sum(NS(1:k))).const<river(k).const
                new_river=stream(j+sum(NS(1:k)));
                stream(j+sum(NS(1:k)))=river(k);
                river(k)=new_river;
                Yes=1;
            end
            
            if Yes==1
                if (river(k).const<=eps && sea.const<=epss && river(k).cost<sea.cost) || (sea.const>epss && river(k).const<=epss)
                    new_sea=river(k);
                    river(k)=sea;
                    sea=new_sea;
                elseif sea.const>epss && river(k).const>epss && sea.const>river(k).const
                    new_sea=river(k);
                    river(k)=sea;
                    sea=new_sea;
                end    
            end          
        end
    end
    % Moving Rivers to Sea
    for j=1:Nsr-1
        river(j).position=river(j).position+2.*rand(1,nvars).*(sea.position-river(j).position);
        
        river(j).position=min(river(j).position,UB);
        river(j).position=max(river(j).position,LB);
        
        river(j).cost=objective_function(river(j).position);
        c=constraints(river(j).position);
        river(j).const=sum(c(c>epss));
        
        if (sea.const>epss && river(j).const<=epss) || (sea.const<=epss && river(j).const<=epss && river(j).cost<sea.cost)
            new_sea=river(j);
            river(j)=sea;
            sea=new_sea;
        elseif sea.const>epss && river(j).const>epss && river(j).const<sea.const
            new_sea=river(j);
            river(j)=sea;
            sea=new_sea;
        end
    end
    % Evaporation condition and raining process
    % Check the evaporation condition for rivers and sea
    for k=1:Nsr-1
        if (norm(sea.position-river(k).position)<dmax || rand<0.1)
            for j=1:NB(k)
                stream(j+sum(NS(1:k))).position=LB+randn(1,nvars).*(UB-LB);
            end
        end
    end
    % Check the evaporation condition for streams and sea
    for j=1:NS(1)
        if norm(sea.position-stream(j).position)<dmax
            stream(j).position=sea.position+sqrt(0.1).*rand(1,nvars);
        end
    end
    
    dmax=dmax-(dmax/max_it);
    
    disp(['Iteration:    ',num2str(i),'    Fmin:    ',num2str(sea.cost),'    Sum_Const:    ',num2str(sea.const)]);
    FF(i)=sea.cost;
end
%% Results and plot
toc;
Elapsed_Time=toc;
plot(FF,'Linewidth',2);
xlabel('Number of Iteration');
ylabel('Function Values');
NFEs=max_it*Npop;
Xmin=sea.position;
Fmin=objective_function(Xmin);
c=constraints(sea.position);
Sum_Const=sum(c(c>epss));

end

