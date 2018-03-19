import boto3

ec2 = boto3.resource('ec2', region_name='us-east-1')

instances = ec2.create_instances(ImageId='ami-d15a75c7',
                                 MinCount=1,
                                 MaxCount=1,
                                 KeyName='eph_cloud',
                                 InstanceType='c4.8xlarge')



