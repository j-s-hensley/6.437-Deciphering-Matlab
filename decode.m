function [] = decode(ciphertext,output_file_name)
%% pre-processing
load(strcat(pwd,'/language_parameters.mat'));
Ml = log(letter_transition_matrix);
Pl = log(letter_probabilities);

%% find the cipher

f = initf(ciphertext,alphabet);
counter = 0;
threshold = 1.5e9/max(length(ciphertext),5e4);
while counter<threshold
    fprime = genf(f);
    a = min(1,comparepfl(fprime,f,ciphertext,alphabet,Ml,Pl));
    u = binornd(1,a);
    if u
        counter = 0;
        f = fprime;
    else
        counter = counter+1;
    end     
end

%% Output

deciphered = invert_text(f,ciphertext,alphabet);
fID = fopen(output_file_name,'w');
fprintf(fID,deciphered);
fclose(fID);
end

%%
function [f] = genf(fin)
%generate a new permutation function based on an input
%according to V = 1/756 iff f' and f differ by two sym assignments, 0 ow
bad = true;
while bad
    randos = ceil(28*rand(1,2));
    if randos(1) ~= randos(2)
       bad = false;
    end
end
a = min(randos);
b = max(randos);
f = [fin(1:(a-1)); fin(b); fin(a+1:b-1); fin(a); fin(b+1:end)];
end
%%
function [expprob] = comparepfl(finnum,finden,y,alphabet,Ml,Pl)

xnumerator = invtext(finnum,y,alphabet);
xdenominator = invtext(finden,y,alphabet);
a = Pl(xnumerator(1)) + sum(Ml(sub2ind([28 28],xnumerator(2:end),xnumerator(1:end-1))));
b = Pl(xdenominator(1)) + sum(Ml(sub2ind([28 28],xdenominator(2:end),xdenominator(1:end-1))));
if a==b
    expprob = 1;
else
    expprob = exp(a-b);
end
end
%%
function fstart = initf(y,alphabet)

indices = cell(29,1);
indices{29} = [];
followers = cell(28,1);
candper = 29;
for i = 1:28
    indices{i} = find(y(1:end-1)==alphabet(i));
    followers{i} = y(indices{i}+1);
    if all(followers{i}==followers{i}(1))&&length(indices{candper})<length(indices{i})
        candper = i;
    end
end
candspace = find(alphabet==followers{candper}(1));
fstart = [1:28]';
fstart(fstart==candspace) = [];
fstart(fstart==candper) = [];
fstart = [fstart; candspace; candper];
end

%%
function [plain] = invert_text(fin,y,alphabet)

finversion = zeros(28,1);
sortedtext = cell(28,1);
for i=1:28
    finversion(i) = find(fin==i);
    sortedtext{i} = find(y==alphabet(i));
end
invalph = alphabet(finversion);
plain = y;
for i=1:28
    plain(sortedtext{i}) = invalph(i);
end
end

%%
function [xnum] = invtext(fin,y,alphabet)

xnum = zeros(length(y),1);
for i=1:28
    xnum(y==alphabet(i)) = find(fin==i);
end
end