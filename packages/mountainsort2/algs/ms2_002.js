#!/usr/bin/env nodejs

var common=require('./common.js');
var os=require('os');

exports.spec=function() {
	var spec0={};
	spec0.name='mountainsort.ms2_002';
	spec0.version='0.1';

	spec0.inputs=[
        {name:"timeseries",description:"preprocessed timeseries (M x N)",optional:false},
        {name:"prescribed_event_times",description:"Timestamps for all events",optional:true},
        {name:"event_times",description:"Timestamps for all events",optional:true},
        {name:"amplitudes",description:"Amplitudes for all events",optional:true},
        {name:"clips",description:"Event clips (from preprocessed timeseries)",optional:true}
    ];
    spec0.outputs=[
    	{name:"event_times_out",optional:true},
    	{name:"amplitudes_out",optional:true},
    	{name:"clips_out",optional:true},
      {name:"firings_out",description:"The labeled events (R x L), R=3 or 4",optional:true}
    ];

	spec0.parameters=[];
	spec0.parameters.push({name:"samplerate",description:"sample rate for timeseries",optional:false});
	spec0.parameters.push({name:'segment_duration_sec',optional:true,default_value:3600});
	spec0.parameters.push({name:'num_threads',optional:true,default_value:0});
	spec0.parameters.push({name:"central_channel",optional:true,default_value:0});
	spec0.parameters.push({name:"clip_size_msec",optional:true,default_value:2});
	spec0.parameters.push({name:"detect_interval_msec",optional:true,default_value:1});
	spec0.parameters.push({name:"detect_threshold",optional:true,default_value:3.0});
	spec0.parameters.push({name:"detect_sign",optional:true,default_value:0});
	spec0.parameters.push({name:"consolidate_clusters",optional:true,default_value:'false'});
	spec0.parameters.push({name:"consolidation_factor",optional:true,default_value:0.9});
	spec0.parameters.push({name:"fit_stage",optional:true,default_value:'false'});
	spec0.parameters.push({name:'subsample_factor',optional:true,default_value:1});
	spec0.parameters.push({name:'channels',optional:true,default_value:''});
	return common.clone(spec0);
};

exports.run=function(opts,callback) {
	var tmpfiles=[]; //all the temporary files to get removed at the end
	opts.temp_prefix=opts.temp_prefix||'00'; //in case the user has specified the temp_prefix
	if (!opts.num_threads) opts.num_threads=os.cpus().length;
	if (!opts.central_channel) opts.central_channel=0;
	var central_channel2=opts.central_channel;
	if ((opts.central_channel)&&(opts.channels)) {
		var channels_list=opts.channels.split(',');
		for (var i=0; i<channels_list.length; i++) {
			if (channels_list[i]==opts.central_channel) {
				central_channel2=(i+1);
			}
		}
	}
	
	if (opts.clips) {
		if (!opts.event_times) {
			console.error('If you specify clips as input, then you must also provide event_times');
			process.exit(-1);	
		}
	}

	if (opts.firings_out) {
		if ((opts.event_times_out)||(opts.amplitudes_out)||(opts.clips_out)) {
			console.error('If you request firings_out you cannot also request event_times_out, clips_out, or amplitudes_out');
			process.exit(-1);
		}
	}
	
	if (!opts._tempdir) {
		console.error('Mandatory parameter _tempdir is empty'); //we need a _tempdir to store all the temporary files
		process.exit(-1);
	}
	if (!opts.samplerate) {
		console.error('opts.samplerate is zero or empty');
		process.exit(-1);	
	}

	console.log ('Using '+opts.num_threads+' threads');

	var steps=[]; //the processing steps
	var info={}; //the info about the size of the arrays
	var segments=[]; //time segments for breaking up the computation
	var event_times=opts.event_times||mktmp('event_times.mda'); //across all segments
	var amplitudes=opts.amplitudes||mktmp('amplitudes.mda'); //across all segments
	var clips=opts.clips||mktmp('clips.mda'); //across all segments (subsampled collection)
	var labels=mktmp('labels.mda'); //across all segments
	var labels2=mktmp('labels2.mda'); //across all segments, after consolidation (labels to remove marked as 0)
	var firings=mktmp('firings.mda'); //across all segments
	var firings_fit=mktmp('firings_fit.mda'); //after fit stage
	//var cluster_metrics=mktmp('cluster_metrics.json');

	var num_intersegment_threads=1; //set later

	///////////////////////////////////////////////////////////////
	steps.push(function(cb) {
		//get the info about the size of the concatenated dataset
		//and set up the segments
		STEP_read_info_from_input_files(function() {
			cb();
		});		
	});
	///////////////////////////////////////////////////////////////
	steps.push(function(cb) {
		//detect events for each segment (if not provided in input)
		//compute the amplitudes
		STEP_presort(function() {
			cb();
		});
	});
	if (!opts.clips) {
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//combine all the event times
			STEP_combine_event_times(function() {
				cb();
			});
		});
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//combine all the event amplitudes
			STEP_combine_amplitudes(function() {
				cb();
			});
		});
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//extract all the clips (subsampled collection)
			STEP_extract_clips(function() {
				cb();
			});
		});
	}
	
	if (opts.firings_out) {
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//the actual clustering
			STEP_sort_clips(function() {
				cb();
			});
		});
		if (opts.consolidate_clusters=='true') {
			///////////////////////////////////////////////////////////////
			steps.push(function(cb) {
				STEP_consolidate_clusters(function() {
					cb();
				});
			});
		}
		else {
			labels2=labels;
		}
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//create firings file
			STEP_create_firings(function() {
				cb();
			});
		});
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//fit stage, etc
			STEP_postsort(function() {
				cb();
			});
		});
		/*
		///////////////////////////////////////////////////////////////
		steps.push(function(cb) {
			//get the cluster metrics
			STEP_cluster_metrics(function() {
				cb();
			});
		});
		*/
	}
	///////////////////////////////////////////////////////////////
	steps.push(function(cb) {
		//write output files
		STEP_write_output_files(function() {
			cb();
		});
	});
	///////////////////////////////////////////////////////////////
	steps.push(function(cb) {
		//remove the temporary files
		STEP_cleanup(function() {
			cb();
		});
	});
	///////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	//Run all steps
	common.foreach(steps,{label:'<processing steps>'},function(ii,step,cb) {
		console.log ('');
		console.log ('--------------------------- SORTING STEP '+(ii+1)+' of '+steps.length +' -----------');
		var timer=new Date();
		step(function() {
			console.log ('SORTING STEP '+(ii+1)+': Elapsed time (sec): '+get_elapsed_sec(timer));
			cb();
		});
	},function() {
		callback();
	});
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

	function STEP_presort(step_callback) {
		var presort_steps=[];
		for (var iseg=0; iseg<segments.length; iseg++) {
			add_presort_step(iseg);
		}
		function add_presort_step(iseg) {
			var segment=segments[iseg];
			segment.timeseries0=mktmp('timeseries_segment_'+iseg+'.mda'); //will be removed after event times obtained
			segment.event_times0=mktmp('event_times0_segment_'+iseg+'.mda');
			segment.event_times1=mktmp('event_times1_segment_'+iseg+'.mda'); //after applying timestamp offset
			segment.amplitudes0=mktmp('amplitudes_segment_'+iseg+'.mda'); //the amplitudes corresponding to the event times
			var intersegment_steps=[];
			presort_steps.push(function(cb) {
				if (!opts.clips) {
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> pre-sort: Extracting timeseries for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						extract_segment_timeseries(opts.timeseries,segment.timeseries0,segment.t1,segment.t2,opts.channels,function() {
							cb2();
						});
					});
					if (opts.prescribed_event_times) {
						//use has supplied the event times as input
						intersegment_steps.push(function(cb2) {
							console.log ('>>>>>>>>>>>> pre-sort: Using input event times for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
							extract_segment_event_times(opts.prescribed_event_times,segment.event_times0,segment.t1,segment.t2,function() {
								cb2();
							});
						});
					}
					else {
						intersegment_steps.push(function(cb2) {
							console.log ('>>>>>>>>>>>> pre-sort: Detect events for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
							detect_events(segment.timeseries0,segment.event_times0,central_channel2,function() {
								cb2();
							});
						});
					}
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> pre-sort: Compute amplitudes for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						compute_amplitudes(segment.timeseries0,segment.event_times0,central_channel2,segment.amplitudes0,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> pre-sort: Applying timestamp offset for events in segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						apply_timestamp_offset(segment.event_times0,segment.event_times1,segment.t1,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> pre-sort: Removing timeseries data to preserve disk space for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						common.remove_temporary_files([segment.timeseries0],function() {
							segment.timeseries0='';
							cb2();
						});
					});
				}
				//Run all intersegment steps
				common.foreach(intersegment_steps,{num_parallel:1,label:'<pre-sort intersegment steps>'},function(ii,step0,cb0) {
					step0(function() {
						cb0();
					});
				},function() {
					cb();
				});
			});
		}
		//Run all process1 steps
		common.foreach(presort_steps,{num_parallel:opts.num_threads,label:'<process1 steps>'},function(ii,step,cb) {
			console.log ('');
			console.log ('--------------------------- PRE-SORT PROCESS SEGMENT '+(ii+1)+' of '+presort_steps.length +' -----------');
			var timer=new Date();
			step(function() {
				console.log ('Elapsed time for process1 '+(ii+1)+' of '+(presort_steps.length)+' (sec): '+get_elapsed_sec(timer));
				cb();
			});
		},function() {
			step_callback();
		});
	}

	function STEP_postsort(step_callback) {
		num_intersegment_threads=Math.floor(opts.num_threads/segments.length);
		if (num_intersegment_threads<1) num_intersegment_threads=1;

		var postsort_steps=[];
		for (var iseg=0; iseg<segments.length; iseg++) {
			add_postsort_step(iseg);
		}
		function add_postsort_step(iseg) {
			var segment=segments[iseg];
			var intersegment_steps=[];
			postsort_steps.push(function(cb) {
				if (opts.fit_stage=='true') {
					segment.timeseries0=mktmp('timeseries_segment_'+iseg+'.mda'); //will be removed after fit
					segment.firings0=mktmp('firings_segment_'+iseg+'.mda');
					segment.firings_fit0=mktmp('firings_fit0_segment_'+iseg+'.mda');
					segment.firings_fit1=mktmp('firings_fit1_segment_'+iseg+'.mda'); //after timestamp offset applied
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> post-sort----: Extracting timeseries for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						extract_segment_timeseries(opts.timeseries,segment.timeseries0,segment.t1,segment.t2,opts.channels,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> post-sort: Extracting firings for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						extract_segment_firings(firings,segment.firings0,segment.t1,segment.t2,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> post-sort: Fit stage for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						fit_stage(segment.timeseries0,segment.firings0,segment.firings_fit0,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> post-sort: Applying timestamp offsets for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						apply_timestamp_offset(segment.firings_fit0,segment.firings_fit1,segment.t1,function() {
							cb2();
						});
					});
					intersegment_steps.push(function(cb2) {
						console.log ('>>>>>>>>>>>> post-sort: Removing timeseries data to preserve disk space for segment '+(iseg+1)+' ('+segment.t1+','+segment.t2+')...\n');	
						common.remove_temporary_files([segment.timeseries0],function() {
							segment.timeseries0='';
							cb2();
						});
					});
				}
				else {
					segment.firings_fit1=segment.firings0;
				}
				
				//Run all intersegment steps
				common.foreach(intersegment_steps,{num_parallel:1,label:'<post-sort intersegment steps>'},function(ii,step0,cb0) {
					step0(function() {
						cb0();
					});
				},function() {
					cb();
				});
			});
		}
		//Run all process2 steps
		common.foreach(postsort_steps,{num_parallel:opts.num_threads,label:'<process2 steps>'},function(ii,step,cb) {
			console.log ('');
			console.log ('--------------------------- POST-SORT PROCESS SEGMENT '+(ii+1)+' of '+postsort_steps.length +' -----------');
			var timer=new Date();
			step(function() {
				console.log ('Elapsed time for process2 '+(ii+1)+' of '+(postsort_steps.length)+' (sec): '+get_elapsed_sec(timer));
				cb();
			});
		},function() {
			if (opts.fit_stage=='true') {
				var firings_fit_list=[];
				for (var iseg=0; iseg<segments.length; iseg++) {
					firings_fit_list.push(segments[iseg].firings_fit1);
				}
				console.log ('Concatenating firings...');
				common.mp_exec_process('mountainsort.concat_firings',
					{firings_list:firings_fit_list},
					{firings_out:firings_fit},
					{},
					step_callback
				);
			}
			else {
				firings_fit=firings;
				step_callback();
			}
		});
	}

	function apply_timestamp_offset(firings,firings_out,timestamp_offset,callback) {
		common.mp_exec_process('mountainsort.apply_timestamp_offset',
			{firings:firings},
			{firings_out:firings_out},
			{timestamp_offset:timestamp_offset},
			callback
		);
	}

	function STEP_combine_event_times(combine_event_times_callback) {
		console.log ('-------------------- COMBINING EVENT TIMES -------------------');
		var event_times_list=[];
		for (var ii in segments) {
			var segment=segments[ii];
			event_times_list.push(segment.event_times1);
		}
		common.mp_exec_process('mountainsort.concat_event_times',
			{event_times_list:event_times_list},
			{event_times_out:event_times},
			{_request_num_threads:opts.num_threads},
			combine_event_times_callback
		);
	}

	function STEP_combine_amplitudes(combine_amplitudes_callback) {
		console.log ('-------------------- COMBINING AMPLITUDES -------------------');
		var amplitudes_list=[];
		for (var ii in segments) {
			var segment=segments[ii];
			amplitudes_list.push(segment.amplitudes0);
		}
		common.mp_exec_process('mountainsort.concat_event_times', //a hack, use the same processor
			{event_times_list:amplitudes_list},
			{event_times_out:amplitudes},
			{_request_num_threads:opts.num_threads},
			combine_amplitudes_callback
		);
	}

	function STEP_extract_clips(extract_clips_callback) {
		console.log ('-------------------- EXTRACTING CLIPS -------------------');
		var clip_size=Math.ceil(opts.clip_size_msec/1000*opts.samplerate);
		common.mp_exec_process('mountainsort.extract_clips',
			{timeseries:opts.timeseries,event_times:event_times},
			{clips_out:clips},
			{clip_size:clip_size,channels:opts.channels||''},
			extract_clips_callback
		);
	}

	function STEP_sort_clips(sort_clips_callback) {
		console.log ('-------------------- SORTING CLIPS -------------------');
		common.mp_exec_process('mountainsort.sort_clips',
			{clips:clips},
			{labels_out:labels},
			{_request_num_threads:opts.num_threads},
			sort_clips_callback
		);	
	}

	function STEP_consolidate_clusters(consolidate_clusters_callback) {
		console.log ('-------------------- CONSOLIDATING CLUSTERS -------------------');
		common.mp_exec_process('mountainsort.consolidate_clusters',
			{timeseries:opts.timeseries,event_times:event_times,labels:labels},
			{labels_out:labels2},
			{
				central_channel:opts.central_channel,
				consolidation_factor:opts.consolidation_factor,
				_request_num_threads:opts.num_threads
			},
			consolidate_clusters_callback
		);	
	}

	function STEP_create_firings(create_firings_callback) {
		console.log ('-------------------- CREATING FIRINGS -------------------');
		common.mp_exec_process('mountainsort.create_firings',
			{event_times:event_times,labels:labels2,amplitudes:(amplitudes||'')},
			{firings_out:firings},
			{central_channel:opts.central_channel,_request_num_threads:opts.num_threads},
			create_firings_callback
		);	
	}		

	function fit_stage(timeseries,firings,firings_out,callback) {
		common.mp_exec_process('mountainsort.fit_stage',
			{timeseries:timeseries,firings:firings},
			{firings_out:firings_out},
			{
			},
			callback
		);
	}

	function STEP_write_output_files(write_output_files_callback) {
        var src=[];
        var dst=[];
        if (opts.event_times_out) {
        	src.push(event_times); dst.push(opts.event_times_out);
        }
        if (opts.amplitudes_out) {
        	src.push(amplitudes); dst.push(opts.amplitudes_out);
        }
        if (opts.clips_out) {
        	src.push(clips); dst.push(opts.clips_out);
        }
        if (opts.firings_out) {
        	src.push(firings_fit); dst.push(opts.firings_out);
        }

		common.copy_files(src,dst,function() {
			write_output_files_callback();
		});
	}

	/*
	function STEP_cluster_metrics(cluster_metrics_callback) {
		////////////////////////////////////////////////////////
		//cluster metrics
		if (opts.cluster_metrics_out) {
			console.log ('>>>>> Cluster metrics');
			common.mp_exec_process('mountainsort.cluster_metrics',
					{timeseries:opts.timeseries,firings:firings_fit},
					{cluster_metrics_out:cluster_metrics},
					{samplerate:opts.samplerate,_request_num_threads:opts.num_threads},
					function() {
						cluster_metrics_callback();
						//done
					}
			);
		}
		else {
			cluster_metrics_callback();
			//done
		}
	}
	*/

	function STEP_cleanup(callback) {
		common.remove_temporary_files(tmpfiles,callback);
	}

	function detect_events(timeseries,event_times_out,central_channel,callback) {
		var detect_interval=Math.ceil(opts.detect_interval_msec/1000*opts.samplerate);
		common.mp_exec_process('mountainsort.detect_events',
			{timeseries:timeseries},
			{event_times_out:event_times_out},
			{
				central_channel:central_channel,
				detect_threshold:opts.detect_threshold,
				detect_interval:detect_interval,
				sign:opts.detect_sign,
				subsample_factor:opts.subsample_factor||1,
				_request_num_threads:num_intersegment_threads	
			},
			callback
		);
	}

	function extract_segment_event_times(event_times,event_times_out,t1,t2,callback) {
		common.mp_exec_process('mountainsort.extract_segment_firings',
			{firings:event_times},
			{firings_out:event_times_out},
			{t1:t1,t2:t2},
			callback
		);
	}

	function compute_amplitudes(timeseries,event_times,central_channel,amplitudes_out,callback) {
		common.mp_exec_process('mountainsort.compute_amplitudes',
			{timeseries:timeseries,event_times:event_times},
			{amplitudes_out:amplitudes_out},
			{
				central_channel:central_channel,
				_request_num_threads:num_intersegment_threads
			},
			callback
		);	
	}

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	

	function mktmp(name) {
		var temp_prefix=opts._temp_prefix||'00';
		var path=opts._tempdir+'/'+temp_prefix+'-'+name;
		tmpfiles.push(path);
		return path;
	}

	function STEP_read_info_from_input_files(callback) {
		console.log ('Reading info from input files...');
		if (typeof(opts.timeseries)!='object') {
			common.read_mda_header(opts.timeseries,function (header) { // Read the .mda header for the timeseries
				info.M=header.dims[0];
				info.N=header.dims[1];
				if (opts.channels) {
					info.M=opts.channels.split(',').length;
				}
				do_create_segments();
			});
		}
		else {
			get_header_for_concatenation_of_timeseries(opts.timeseries,function(header) {
				info.M=header.dims[0];
				info.N=header.dims[1];
				if (opts.channels) {
					info.M=opts.channels.split(',').length;
				}
				do_create_segments();
			});
		}
		function do_create_segments() {
			var segment_duration=Math.ceil(opts.segment_duration_sec*opts.samplerate);
			segments=create_segments(info.N,segment_duration,segment_duration);
			if (segments.length===0) {
				console.log ('Error: no segments created (N='+info.N+')');
				process.exit(-1);
			}
			num_intersegment_threads=Math.floor(opts.num_threads/segments.length);
			if (num_intersegment_threads<1) num_intersegment_threads=1;
			callback();
		}
	}

	function extract_segment_timeseries(ts_in,ts_out,t1,t2,channels,callback) {
		common.grab_lock('extract_segment_timeseries',function() {
			common.mp_exec_process('mountainsort.extract_segment_timeseries',
				{timeseries:ts_in},
				{timeseries_out:ts_out},
				{t1:t1,t2:t2,channels:(channels||'')},
				function() {
					common.release_lock('extract_segment_timeseries');
					callback();
				}
			);
		});
	}

	function extract_segment_firings(firings_in,firings_out,t1,t2,callback) {
		common.mp_exec_process('mountainsort.extract_segment_firings',
			{firings:firings_in},
			{firings_out:firings_out},
			{t1:t1,t2:t2},
			callback
		);
	}
};

function get_elapsed_sec(timer) {
	var stamp=new Date();
	return (stamp-timer)/1000;
}

function create_segments(N,segment_size,shift_size) {
	var ret=[];
	for (var i=0; i<N; i+=shift_size) {
		var t1=i,t2=t1+segment_size-1;
		if (t2>=N) {
			t2=N-1
			ret.push({t1:t1,t2:t2});
			break;
		}
		else {
			ret.push({t1:t1,t2:t2});
		}
	}
	return ret;
}

function get_header_for_concatenation_of_timeseries(ts_list,callback) {
	if (ts_list.length==0) {
		console.error('ts_list is empty.');
		process.exit(-1);
	}
	common.read_mda_header(ts_list[0],function(header0) {
		header0.dims[1]=0;
		common.foreach(ts_list,{},function(ii,ts,cb) {
			common.read_mda_header(ts,function(header1) {
				header0.dims[1]+=header1.dims[1];
				cb();
			});
		},function() {
			callback(header0);
		});
	});
	
}
