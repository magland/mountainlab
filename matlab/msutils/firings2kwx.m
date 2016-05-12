function firings2kwx(timeseries,firings,outfile,T,num_features_per_channel)

% Notes
% set masks to all 1's for now

if nargin<1, test_firings2kwx; return; end;

if nargin<4
    T=41;
end;
if nargin<5
    num_features_per_channel=3;
end;


mfile_path=fileparts(mfilename('fullpath'));
template=[mfile_path,'/template.kwik'];

times=firings(2,:);
labels=firings(3,:);

M=size(timeseries,1);
N=size(timeseries,2);
K=max(labels);
L=length(times);

if (exist(outfile,'file')) delete(outfile); end;

num_features=num_features_per_channel*M;
features_masks=zeros(2,num_features,L);
features_masks(2,:)=1; %set masks to all 1's for now
for m=1:M
    clips=ms_extract_clips2(timeseries(m,:),times,T,1);
    features_masks(1,(m-1)*num_features_per_channel+1:m*num_features_per_channel,:)=ms_event_features(clips,num_features_per_channel);
end;

clips=ms_extract_clips2(timeseries,times,T,1);

%scale to int16
maxabs=max(abs(clips(:)));
clips=clips/maxabs*(32767-1);

h5create(outfile,'/channel_groups/0/features_masks',size(features_masks),'Datatype','single');
h5write(outfile,'/channel_groups/0/features_masks',features_masks);
h5create(outfile,'/channel_groups/0/waveforms_filtered',size(clips),'Datatype','int16','ChunkSize',[M,T,min(195,L)]);
h5write(outfile,'/channel_groups/0/waveforms_filtered',clips);
h5create(outfile,'/channel_groups/0/waveforms_raw',size(clips),'Datatype','int16','ChunkSize',[M,T,min(195,L)]);
h5write(outfile,'/channel_groups/0/waveforms_raw',clips);

h5writeatt(outfile,'/','TITLE','');
h5writeatt(outfile,'/','CLASS','GROUP');
h5writeatt(outfile,'/','VERSION','1.0');
h5writeatt(outfile,'/','PYTABLES_FORMAT_VERSION','2.1');

h5writeatt(outfile,'/channel_groups','TITLE','');
h5writeatt(outfile,'/channel_groups','CLASS','GROUP');
h5writeatt(outfile,'/channel_groups','VERSION','1.0');

end

function test_firings2kwx

M=2;
N=100;
K=2;
L=20;
times=1:L;
labels=randi(K,size(times));
firings=zeros(3,L); firings(2,:)=times; firings(3,:)=labels;
timeseries=randn(M,N);

firings2kwx(timeseries,firings,'test1.kwx');
h5disp('test1.kwx');


end