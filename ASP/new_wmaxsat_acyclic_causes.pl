%%% Modified version of Hyttinen's algorithm: we just add the causal statements.'

causes(X,Y) :- edge(X,Y), X!=Y, node(X), node(Y).
causes(X,Y) :- causes(X,Z), causes(Z,Y), X!=Y, X!=Z, Y!=Z, node(X), node(Y), node(Z).
:- causes(X,X), node(X).
:- causes(X,Y), causes(Y,X), node(X), node(Y), X!=Y.

%% { indep(X,Y,Z) } :- node(X), node(Y), set(Z), not ismember(Z, X), not ismember(Z,Y), X != Y.
%% :- not indep(Y,X,Z), indep(X,Y,Z), node(X), node(Y), set(Z), not ismember(Z, X), not ismember(Z,Y), X != Y.
%% :- th(X,Y,C,J,M), indep(X,Y,C), X != Y.
%% :- th(Y,X,C,J,M), indep(X,Y,C), X != Y.
%% :- hh(X,Y,C,J,M), indep(X,Y,C), X != Y.
%% :- tt(X,Y,C,J,M), indep(X,Y,C), X != Y.
%% :- not th(X,Y,C,J,M), not th(Y,X,C,J,M), not hh(X,Y,C,J,M), not tt(X,Y,C,J,M), not indep(X,Y,C), X != Y.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{ before(X,Y) } :- node(X), node(Y), X != Y.
{ edge(X,Y) } :- node(X), node(Y), before(X,Y), X != Y.
{ conf(X,Y) } :- node(X), node(Y), X<Y.

%the order is  transitive
:-not before(X,Z), before(X,Y), before(Y,Z), node(X), node(Y), node(Z), X != Y, Y != Z, X != Z.
%before relation cannot happen in both directions
:-not before(Z,Y), not before(Y,Z), node(Y), node(Z), Y != Z.
:-before(Z,Y), before(Y,Z), node(Y), node(Z), Y != Z.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOTTOM

th(X,Y,0,0,0) :- edge(X,Y).
hh(X,Y,0,0,0) :- conf(X,Y).

%notice that acyclicity can be done only when cset=0
%not needed for anything!
%:- tt(X,Y,0,J,M). %no tail tail paths if no conditioning

%requirements of intervention operation on z
%is_member($z,$J) && !is_member($z,$C) && !is_member($z,$M) )


%%%%%%%%%%%INTERVENTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%notice that here must allow tt(X,X,...) paths
tt(X,Y,C,J,M):- tt(X,Y,C,Jsub,M), 
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				intervene(C,Jsub,Z,J,M).

th(X,Y,C,J,M):- th(X,Y,C,Jsub,M), 
				X != Y,
				not ismember(J,Y),
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				intervene(C,Jsub,Z,J,M).

hh(X,Y,C,J,M):- hh(X,Y,C,Jsub,M), 
				X < Y,
				not ismember(J,Y), not ismember(J,X),
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				intervene(C,Jsub,Z,J,M).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%CONDITION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt(X,Y,C,J,M):- tt(X,Y,Csub,J,M), 
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).

tt(X,Y,C,J,M):- th(X,Z,Csub,J,M),th(Y,Z,Csub,J,M), 
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th(X,Y,C,J,M):- th(X,Y,Csub,J,M), 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).

th(X,Y,C,J,M):- th(X,Z,Csub,J,M),
				{ hh(Z,Y,Csub,J,M); hh(Y,Z,Csub,J,M) } >= 1,
				X != Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).
				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hh(X,Y,C,J,M):- hh(X,Y,Csub,J,M), 
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).

%three options for symmetry ZXY,XZY,XYZ				
hh(X,Y,C,J,M):- { hh(Z,X,Csub,J,M); hh(X,Z,Csub,J,M) } >= 1,
				{ hh(Z,Y,Csub,J,M); hh(Y,Z,Csub,J,M) } >= 1,
				X < Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				condition(Csub,Z,C,J,M).
				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%MARGINALIZATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt(X,Y,C,J,M):- tt(X,Y,C,J,Msub), 
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

tt(X,Y,C,J,M):- th(X,Z,C,J,Msub), 
				{ tt(Z,Y,C,J,Msub); tt(Y,Z,C,J,Msub) } >= 1,
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

tt(X,Y,C,J,M):- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1, 
				th(Y,Z,C,J,Msub),
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

tt(X,Y,C,J,M):- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1, 
				{ tt(Z,Y,C,J,Msub); tt(Y,Z,C,J,Msub) } >= 1,
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

tt(X,Y,C,J,M):- th(X,Z,C,J,Msub),th(Y,Z,C,J,Msub),tt(Z,Z,C,J,Msub), 
				X <= Y, 
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th(X,Y,C,J,M):- th(X,Y,C,J,Msub), 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

th(X,Y,C,J,M):- th(X,Z,C,J,Msub),th(Z,Y,C,J,Msub), 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

th(X,Y,C,J,M):- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1,
				{ hh(Z,Y,C,J,Msub); hh(Y,Z,C,J,Msub) } >= 1, 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

th(X,Y,C,J,M):- { tt(X,Z,C,J,Msub); tt(Z,X,C,J,Msub) } >= 1,
				th(Z,Y,C,J,Msub), 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

th(X,Y,C,J,M):- th(X,Z,C,J,Msub), 
				tt(Z,Z,C,J,Msub),
				{ hh(Z,Y,C,J,Msub); hh(Y,Z,C,J,Msub) } >= 1, 
				X != Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hh(X,Y,C,J,M):- hh(X,Y,C,J,Msub),
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).


hh(X,Y,C,J,M):- { hh(X,Z,C,J,Msub); hh(Z,X,C,J,Msub) } >= 1,
				th(Z,Y,C,J,Msub),
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

hh(X,Y,C,J,M):- th(Z,X,C,J,Msub),
				th(Z,Y,C,J,Msub),
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),				
				marginalize(C,J,Msub,Z,M).

hh(X,Y,C,J,M):- th(Z,X,C,J,Msub),
				{ hh(Y,Z,C,J,Msub); hh(Z,Y,C,J,Msub) } >= 1,
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).

hh(X,Y,C,J,M):- { hh(X,Z,C,J,Msub); hh(Z,X,C,J,Msub) } >= 1,
				{ hh(Y,Z,C,J,Msub); hh(Z,Y,C,J,Msub) } >= 1,
				tt(Z,Z,C,J,Msub),
				X < Y,
				not ismember(C,X), not ismember(C,Y),
				not ismember(M,X), not ismember(M,Y),
				node(X),node(Y),
				marginalize(C,J,Msub,Z,M).
				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fail(X,Y,C,J,M,W) :- th(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- th(Y,X,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- hh(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.
fail(X,Y,C,J,M,W) :- tt(X,Y,C,J,M), indep(X,Y,C,J,M,W), X<Y.

fail(X,Y,C,J,M,W) :- not th(X,Y,C,J,M), not th(Y,X,C,J,M), not hh(X,Y,C,J,M), not tt(X,Y,C,J,M), dep(X,Y,C,J,M,W), X<Y.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#minimize{W,X,Y,C,J,M:fail(X,Y,C,J,M,W) }.

#show .
#show causes/2.
#show edge/2.
#show conf/2.
#show node/1.
