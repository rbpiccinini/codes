#!/bin/sh
host=ec2-46-137-140-90.eu-west-1.compute.amazonaws.com
scp -i ~/master/amazon/piccinini.pem ubuntu@$host:/home/ubuntu/master/sydney/case/VTK.tar.gz .
tar xvfz VTK.tar.gz
