pipeline.py /data/ckilpatrick/spitzer/2017gfo/ch2 197.4487465 -23.3839793 --band ch2 --object AT2017gfo --mopex-dir /data/software/mopex --init-date 2017-08-17 --max-date 2018-01-01 --sn-offset 0.001591598 0.0024214555 --stamp-size 29 > /data/ckilpatrick/spitzer/2017gfo/ch2/pipeline.out &
pipeline.py /data/ckilpatrick/spitzer/2017gfo/ch1 197.4487465 -23.3839793 --band ch1 --object AT2017gfo --mopex-dir /data/software/mopex --init-date 2017-08-17 --max-date 2018-01-01 --sn-offset 0.001591598 0.0024214555 --stamp-size 29 > /data/ckilpatrick/spitzer/2017gfo/ch1/pipeline.out &

pipeline.py /home/cdk8313/project/spitzer/2017gfo/ch1 197.4487465 -23.3839793 --band ch1 --object AT2017gfo --mopex-dir /data/software/mopex --init-date 2017-08-17 --max-date 2018-01-01 --sn-offset 0.001591598 0.0024214555 --stamp-size 29 --nprocesses 120 > /home/cdk8313/project/spitzer/2017gfo/ch1/pipeline.out &

pipeline.py /data/ckilpatrick/spitzer/2017gfo/test 197.4487465 -23.3839793 --download /data/ckilpatrick/spitzer/2017gfo/raw --band ch2 --object AT2017gfo --mopex-dir /data/software/mopex --init-date 2017-08-17 --max-date 2018-01-01 --sn-offset 0.001591598 0.0024214555 --stamp-size 29
