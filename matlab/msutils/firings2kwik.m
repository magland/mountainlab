function firings2kwik(firings,outfile,channel_positions,samplerate)

% Notes
% copying over spikedetekt application data
% adjacency_graph is empty
% using 'Good' cluster group for every cluster
% all same color (3)
% time_fractional all set to zero
% adjusting to zero-based labels? (how about times)

if nargin<1, test_firings2kwik; return; end;

mfile_path=fileparts(mfilename('fullpath'));
template=[mfile_path,'/template.kwik'];

if (ischar(firings))
    firings=readmda(firings);
end;

times=firings(2,:);
labels=firings(3,:);

if (exist(outfile,'file'))
    delete(outfile);
end;

M=size(channel_positions,1);
K=max(labels);
L=length(times);

adjacency_graph=zeros(2,0);
channel_names=cell(1,M);
for ch=1:M
    channel_names{ch}=sprintf('%d',ch-1);
end;
dat_path=outfile;

grp1='/application_data/spikedetekt';
grp2='/application_data/spikedetekt';
h5copy(template,outfile,grp1,grp2);

for ch=1:M
    grp1=sprintf('/channel_groups/0/channels/0');
    grp2=sprintf('/channel_groups/0/channels/%d',ch-1);
    h5copy(template,outfile,grp1,grp2);
    h5writeatt(outfile,grp2,'name',channel_names{ch});
    h5writeatt(outfile,grp2,'position',channel_positions(ch,:));
end;
h5writeatt(outfile,'/channel_groups/0','channel_order',0:M-1);
h5writeatt(outfile,'/channel_groups/0','adjacency_graph',adjacency_graph);

grp1=sprintf('/channel_groups/0/cluster_groups/main');
grp2=sprintf('/channel_groups/0/cluster_groups/main');
h5copy(template,outfile,grp1,grp2);

grp1=sprintf('/channel_groups/0/cluster_groups/original');
grp2=sprintf('/channel_groups/0/cluster_groups/original');
h5copy(template,outfile,grp1,grp2);

for k=1:K
    grp1=sprintf('/channel_groups/0/clusters/main/0');
    grp2=sprintf('/channel_groups/0/clusters/main/%d',k-1);
    h5copy(template,outfile,grp1,grp2);
    h5writeatt(outfile,grp2,'cluster_group',int16(2));
    h5writeatt(outfile,grp2,'color',int16(3));
end;

h5create(outfile,'/channel_groups/0/spikes/recording',L,'Datatype','int32');
h5write(outfile,'/channel_groups/0/spikes/recording',zeros(L,1));
h5create(outfile,'/channel_groups/0/spikes/time_fractional',L,'Datatype','uint8');
h5write(outfile,'/channel_groups/0/spikes/time_fractional',zeros(L,1));
h5create(outfile,'/channel_groups/0/spikes/time_samples',L,'Datatype','uint64');
h5write(outfile,'/channel_groups/0/spikes/time_samples',times);
h5create(outfile,'/channel_groups/0/spikes/clusters/main',L,'Datatype','int32');
h5write(outfile,'/channel_groups/0/spikes/clusters/main',labels);
h5create(outfile,'/channel_groups/0/spikes/clusters/original',L,'Datatype','int32');
h5write(outfile,'/channel_groups/0/spikes/clusters/original',labels);
%h5write(outfile,'/channel_groups/0/spikes/time_fractional',zeros(L,1));
%h5write(outfile,'/channel_groups/0/spikes/time_samples',times0);
%h5write(outfile,'/channel_groups/0/spikes/clusters/main',labels0-1);
%h5write(outfile,'/channel_groups/0/spikes/clusters/original',labels0-1);


grp1=sprintf('/recordings/0');
grp2=sprintf('/recordings/0');
h5copy(template,outfile,grp1,grp2);
h5writeatt(outfile,grp2,'sample_rate',samplerate);
h5writeatt(outfile,[grp2,'/raw'],dat_path,dat_path);

end

function test_firings2kwik

M=2;
K=2;
L=10;
times=1:L;
labels=randi(K,size(times));
positions=rand(M,2);
samplerate=30000;
firings=zeros(3,L); firings(2,:)=times; firings(3,:)=labels;

firings2kwik(firings,'test1.kwik',positions,samplerate);
h5disp('test1.kwik');


end