# Rosetta
Rosetta 文件夹下包含两个子文件夹：cmdline/, code/。

## cmdline/
该文件夹中包含Rosetta Simple Cycpep的使用手册，帮助用户在linux环境下通过command line的形式进行单个序列预测。内容涉及参数说明、单/多进程运行等。

## code/
该文件夹中包含两个py文件：preparation.py, analyze.py。

### preparation.py

集成了pdb文件(pdbid.chain.pdb)/fasta文件批量转换txt文件、预测脚本生成、提取脚本生成。\
使用方法： python preparation.py `pdb文件夹/fasta文件名` `txt文件夹` `输出silent文件夹` `预测脚本文件名` `提取脚本文件名`\
sample: 
```console
python preparation.py /home/ictsun/pdbfiles /home/ictsun/txtfiles /home/ictsun/output_silents /home/ictsun/predict_shell.sh /home/ictsun/extract_shell.sh
```
**注：** 路径最好使用绝对路径，这样可以在任意文件夹下启动`predict_shell.sh`。

获得`predict_shell.sh`之后，通过`sh predict_shell.sh`运行脚本，可以在先前设置的对应的文件夹中看到结果。\
得到silent文件之后，**请新建一个文件夹用于存放从silent中提取的PDB文件**，并且进入该PDB文件夹之后，再运行`sh extract_shell.sh`，可以得到每一个序列的最低能量结构文件。

### analyze.py
集成了交叉比对结构、结果输出为csv、可视化csv。\
交叉比对包括序列对各自对比预测结构和原结构，以及互相比对原/预测结构。\
使用方法： python analyze.py `原pdb文件夹` `序列对id txt文件名` `预测pdb文件夹` `输出结果csv文件`\
其中`序列对id txt文件名`格式为一行一对pdbid与chain，如：
```editorconfig
1a1p.A, 1a0m.A
2ilp.A, 5xco.B
```