#一、最原始的安装方式
install.packages("包名")
library(包名)

#二、利用BiocManager管理R包
#Bioconductor摒弃了此前biocLite的安装方法，提供了BiocManager管理R包的方法，前者经常断点且不能续传，后者显然更稳定一些：
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  #判断BiocManager环境及安装
BiocManager::install("包名")

#三、从本地安装R包
#科学上网的问题一直困扰着科研工作者，如果上述的安装方式均无法成功执行，可以尝试“本地安装”的方式。
#1、首先，各位可以通过以下网址搜索并下载R包源码文件，貌似这个界面并没有检索模块，自行control + F检索即可。
#https://mirrors.tuna.tsinghua.edu.cn/CRAN/
#2、package的源码下载完毕后可在Rstudio工具栏中使用：
#tools→install.packages→browser进行安装：

#在服务器中无GUI的同学可以使用以下命令进行安装：
#install.packages("路径/源码文件", repos = NULL)
#注：以上方法均可以批量安装R包
#把packages名称以c("packages1","packages2","packages3"...)的形式列出即可

#四、ubuntu及服务器下的安装方式
#由于服务器下的权限问题，以上的方法在linux下的R语言中可能出现permission deny的情况，
#因此使用conda 管理R包显然是更合理的选择。具体操作方法如下：
#wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-5.3.1-Linux-x86_64.
#shbash Anaconda3-5.3.1-Linux-x86_64.sh 
#下载并安装condaconda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/conda config 
#--add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/conda config 
#--add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/conda config 
#--add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
#添加conda国内镜像conda search 包名conda install -y bioconductor-包名
#安装之前的检索可以确保R包的名称没有输入错误，毕竟R包的“学名”可能与我们所默认的有所出入