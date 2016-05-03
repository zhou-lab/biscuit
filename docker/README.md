This explains how to run the dockerized version of biscuit. 

For the biscuit project repo see http://github.com/zwdzwd/biscuit
or the biscuit documentation/wiki: http://github.com/zwdzwd/biscuit/wiki

This docker project only aims to create a dockerized container for biscuit. For anything to do with biscuit itself, please see the biscuit project or contact the biscuit developers.
 

Using biscuit with docker

1.) get the container 

	docker pull zackramjan/biscuit

2.) prepare your data
since the container can only see explicitely specified dirs, put your data and reference in a directoy that will be passed to the docker container
for example, I create a "mydata" dir and put my fastq's (R1.fastq and R2.fastq) and ref genome (NC_001416.fa) in it.

	[user@vai ~]# ls mydata
	NC_001416.fa r1.fastq r2.fastq

3.) start the container.  you will need to specifiy the local data dir(for examlple "mydata") with your refs and fastq

	docker run -ti -d --name=biscuit_docker -v mydata:/data -w /data zackramjan/biscuit

4.) use biscuit. you can now perform alignemnts etc, remember that biscuit command are executed as if you're within your specified "data" directory (ie "mydata" in our example)
for example: first I index my genome

	docker exec -it biscuit_docker biscuit index NC_001416.fa

then I perform an alignment. The resulting sam is printed to stdout

	docker exec -it biscuit_docker biscuit align NC_001416.fa r1.fastq r2.fastq 

samtools is also installed within the container, and we can use it to create a bam. Since we are specifying a pipe (which is a shell operation) we enclose it in quotes and pass it bash: 

	docker exec -it biscuit_docker bash -c "biscuit align NC_001416.fa r1.fastq r2.fastq | samtools view -b -o myout.bam"

the result is "myout.bam" which will get written to the data dir you specified ("mydata" in this example")

5.) optionally, you can open an interactive shell from within the container and do your work there.

	docker exec -it biscuit_docker /bin/bash

/data will contain your data directory.
