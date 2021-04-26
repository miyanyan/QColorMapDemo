#ifndef PLOT_H
#define PLOT_H

#include <QWidget>
#include <QTimer>
#include <QHash>
#include "qcustomplot.h"
#include <QDebug>

class HotPlot : public QWidget
{
    Q_OBJECT

public:
    HotPlot(QWidget *parent = nullptr, int contourCnt = 100);
    ~HotPlot();
    //设置要显示的通道
    void setEnabledChannels(QStringList& channels);
    void setEnabledChannels(std::vector<int>& channels);
    //更新每个通道的数据
    void setData(std::vector<double>& channels);
    //设置颜色的显示范围
    void setColorRange(double lower, double upper);
    //是否显示通道名
    void showChannelName(bool flag);
    //是否显示通道圆
    void showChannelCircle(bool flag);
private:
    void initWindow();
    void readChannelAxisFile(QString path);
    void drawElements();
    void test();
    void handleTimeout();
    //根据每个点的位置计算矩阵值
    void calMatrix(std::vector<double> &channels);//DoubleLinear方式
    double doubleLinear(int m, int n, int X, int Y, double V1, double V2, double V3, double V4);
    double singleLinear(int m, int X, double V1, double V2);
    //提前计算
    void precompute();

    //绘图
    QCustomPlot *m_customPlot;
    QCPColorMap *m_colorMap;
    QCPColorScale *m_colorScale;
    //绘制头、鼻、通道
    std::vector<QPointF> m_channelAxis;//64个通道
    QHash<QString, int> m_channelIndex;//64个通道
    QHash<QString, QPointF> m_noseAxis;//鼻尖
    std::vector<QCPItemEllipse*> m_channelCircle;//通道圆
    std::vector<QCPItemText*> m_channelName;//通道名称
    //实际显示的通道索引
    std::vector<int> m_usedChannels;

    //预计算
    std::vector<std::vector<std::vector<double>>> m_channelWeight;//64个通道在每个坐标点的权重系数 [nx][ny][channel]
    std::vector<std::vector<bool>> m_inCircle;//nx * ny 矩阵 的点是否在圆内

    std::vector<std::vector<double>> m_matrix;
    std::vector<std::vector<double>> m_matrix1;

    int m_contourCnt;//轮廓线的个数

    QTimer m_timer;
};
#endif // PLOT_H
