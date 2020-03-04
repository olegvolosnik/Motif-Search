#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;


//FUNCTIONS

void window_quality();
void seq_reader();
void seq_show();
void quality_reader();
void quality_show();
void vertex_generator();
void edge_generator();
void clique_family();
void print_seq_uni_code();
void print_seq_name_or_uni_code();
void most_common();
void clique_family();




//Global variables:
int window, quality_lvl, deletion_nr, data_set_nr;
int tab_uni_code[300][3];


vector <vector <char> > seq;
vector <vector <int> > quality;

class Vertex{
public:

    char v_name[7];
    int v_qual[7];
    int v_id, v_position, v_del_counter;
    int v_uni_code = 0;

    Vertex(const char seq_name[], const int seq_qual[], int seq_id, int seq_position, int deletion_counter){
        for(int i=0; i<window; i++){
            this->v_name[i] = seq_name[i];
            this->v_qual[i] = seq_qual[i];
        }
        this->v_id = seq_id;
        this->v_position = seq_position;
        this->v_del_counter = deletion_counter;
    };
    void show_vertex(){
        for(int i=0; i<window; i++){
            cout<<this->v_name[i];
        }
        for(int i=0; i<window; i++){
            cout<<' '<<this->v_qual[i]<<' ';
        }
        cout<<"seq id = "<<this->v_id<<", seq pos = "<<this->v_position<<", deletion counter = "<<this->v_del_counter<<", uni code = "<<this->v_uni_code<<"  "<<endl;
    }
};

vector<vector<Vertex>> vertex_matrix;

void window_quality(){
    cout<<"\nPlease set a window in range 4-7: ";
    while(true){
        cin>>window;
        if(window > 7 || window < 4){
            cout<<"Please set proper value for window: ";
        }else break;
    }
    if(window>5){
        deletion_nr = 2;
    }else{
        deletion_nr = 1;
    }
    cout<<"\nPlease set a quality in range 1-30: ";
    while(true){
        cin>>quality_lvl;
        if(quality_lvl > 30 || quality_lvl < 1){
            cout<<"Please set proper value for window: ";
        }else break;
    }
    cout<<"window is equal: "<<window<<", and deltion_nr is: "<<deletion_nr<<", and quality is:"<<quality_lvl<<endl;


}
void seq_reader(){
    char header = '>';
    string line;
    int counter = -1;
    fstream file;
    string data_set = "fasta" + to_string(data_set_nr) + ".txt";
    cout<<data_set<<endl;
    file.open(data_set, ios::in);//change on fasta later
    while (!file.eof()){
        getline(file, line);
        if(line[0] == header){
            seq.emplace_back(vector<char>());
            counter++;
            continue;
        }else{
            for(auto c : line){
                if(c != ' ') {
                    seq[counter].push_back(c);
                }
            }
        }
    }
    file.close();
}
void seq_show(){
    cout<<"Sequences: "<<endl;
    for(int i = 0; i < seq.size(); i++){
        cout<<"\nseq["<<i<<"].size() == "<<seq[i].size()<<endl;
        for (char j : seq[i])        {
            cout << j << " ";
        }
    }
}
void quality_reader(){
    char header = '>';
    string line;
    int counter = -1;
    fstream file;
    string data_set = "qual_" + to_string(data_set_nr) + ".txt";
    cout<<data_set<<endl;
    file.open(data_set, ios::in);
    while(!file.eof()){
        getline(file, line);
        if(line[0] == header){
            quality.emplace_back(vector<int>());
            counter++;
            continue;
        }else{
            stringstream ss(line);
            int number;
            while (ss >> number){
                quality[counter].push_back(number);
            }
        }
    }
}
void quality_show(){
    cout<<"\nQualities: "<<endl;
    for(int i = 0; i < quality.size(); i++) {
        cout << "seq[" << i << "]size() == " << quality[i].size() << endl;
        for (int j : quality[i]) {
            cout << j << " ";
        }
        cout << endl;
    }
}
void print_seq_uni_code(){
    //    PRINT SEQ'S OF UNI_CODE
    for (int x = 0; x < seq.size(); x++) {
        for (int y = 0; y < (seq[x].size() - window + 1); y++){
            if(vertex_matrix[x][y].v_uni_code>99)
                cout<<" "<<vertex_matrix[x][y].v_uni_code;
            else if((vertex_matrix[x][y].v_uni_code>9)&&(vertex_matrix[x][y].v_uni_code<100))
                cout<<"  "<<vertex_matrix[x][y].v_uni_code;
            else if(vertex_matrix[x][y].v_uni_code<=9)
                cout<<"   "<<vertex_matrix[x][y].v_uni_code;
        }
        cout<<endl<<endl;
    }
}
void print_seq_name_or_uni_code(){
        for (int x = 0; x < seq.size(); x++) {
        for (int y = 0; y < (seq[x].size() - window + 1); y++){
            if(vertex_matrix[x][y].v_uni_code == 0){
                for(int i = 0; i<window; i++)
                    cout<<vertex_matrix[x][y].v_name[i];
                cout<<" ";
            }else{
                cout<<vertex_matrix[x][y].v_uni_code<<" ";
            }

        }
        cout<<endl;
    }
}


void vertex_generator(){
    int subqual_tab[window];
    char subseq_tab[window];
    int deletion_counter;


    for(int i = 0; i < seq.size(); i++) {
        vertex_matrix.emplace_back();

        for (int k = 0; k < (seq[i].size() - window + 1); k++) {
            deletion_counter = 0;

            for (int m = 0; m < window; m++) {
                subqual_tab[m] = quality[i][k + m];
                if (subqual_tab[m] >= quality_lvl) {
                    subseq_tab[m] = seq[i][k + m];
                } else {
                    subseq_tab[m] = '_';
                    deletion_counter++;
                }
            }
            int seq_id = i, seq_position = k;
            Vertex new_vertex(subseq_tab, subqual_tab, seq_id, seq_position, deletion_counter); //GENERATING NEW OBJECTS
            vertex_matrix[i].push_back(new_vertex);
//            vertex_matrix[i][k].show_vertex();

//            SEQ MOVEMENT PRINTING
//            cout<<"\nseq id = "<<seq_id<<", seq pos = "<<seq_position<<"  ";
//            for(int l = 0; l<=k; l++){if(l>0) cout<<" ";}
//
//            for(int m=0; m<window; m++){cout<<subseq_tab[m];}
//            cout<<" ";
//            for(int m=0; m<window; m++){cout<<subqual_tab[m]<<"  ";}
//            cout<<endl;

        }
    }
}


void edge_generator() {
    bool used_uni_code = false;
    bool seq_eq = false;
    int uni_code = 1;


    for (int x = 0; x < seq.size(); x++) {
        for (int y = 0; y < (seq[x].size() - window + 1); y++) {



            if( vertex_matrix[x][y].v_del_counter > 1)
                continue;  //IF WE HAVE MORE DELETION

            if( vertex_matrix[x][y].v_uni_code != 0)
                continue;  //IF NODE IS ALREADY CHECKED BY NODE-BROTHERS TURN

            if(used_uni_code){
                uni_code++;    used_uni_code = false;}       //GENERATING NEW FREE UNI_CODE



            for (int tmp_x = 0; tmp_x < seq.size(); tmp_x++) {
                for (int tmp_y = 0; tmp_y < (seq[tmp_x].size() - window + 1); tmp_y++){

                    if((x == tmp_x)&&(y == tmp_y))
                        continue;  //IF COMPARE ONE NODE
                    if( vertex_matrix[tmp_x][tmp_y].v_del_counter > deletion_nr)
                        continue;  //IF WE HAVE MORE DELETION



                    seq_eq = true;
                    for(int i = 0; i < window; i++){

                        if(vertex_matrix[x][y].v_name[i] == vertex_matrix[tmp_x][tmp_y].v_name[i])
                            continue;
                        else if(vertex_matrix[tmp_x][tmp_y].v_name[i] == '_')
                            continue;
                        else if(vertex_matrix[x][y].v_name[i] == '_')
                            continue;
                        else{   seq_eq = false;     break;}
                    }


                    if(seq_eq){
                        vertex_matrix[x][y].v_uni_code = uni_code;
                        vertex_matrix[tmp_x][tmp_y].v_uni_code = uni_code;

                        tab_uni_code[uni_code][0]++;
                        tab_uni_code[uni_code][1] = x;
                        tab_uni_code[uni_code][2] = y;
                        used_uni_code = true;

                    }

                }
            }
        }
    }


    most_common();
}


void most_common(){
    int used[seq.size()];
    for(int & i : used)
        i = 0;

    bool seq_eq = false;
    int biggest = 0;
    for(int i = 0; i < 300; i++){
        if(tab_uni_code[i][0] >tab_uni_code[biggest][0])
            biggest = i;
    }

    int x1, y1;
    cout<<"\nBiggest unicode = "<<biggest<<", and it used in matrix "<<tab_uni_code[biggest][0] + 1<<" times."<<endl;
    cout<<"\nX pos = "<<tab_uni_code[biggest][1]<<", Y pos = "<<tab_uni_code[biggest][2]<<endl;
    int new_X = tab_uni_code[biggest][1], new_Y = tab_uni_code[biggest][2];

    //Printing first subseq
    cout<<"X = "<<new_X<<", Y = "<<new_Y<<"\t\t";
    for(int i = 0; i < window; i++) {
        cout << vertex_matrix[new_X][new_Y].v_name[i] << " ";
    }                                   cout<<endl;
    used[new_X] = 1;

    //overwriting one more time all nodes with similar name and printing them
    for (int x = 0; x < seq.size(); x++){
        for (int y = 0; y < (seq[x].size() - window + 1); y++){

            if((x == new_X)&&(y == new_Y))                               continue;  //IF COMPARE ONE NODE
            if( vertex_matrix[x][y].v_del_counter > deletion_nr)         continue;  //IF WE HAVE MORE DELETION

            seq_eq = true;
            for(int i = 0; i < window; i++){

                if(vertex_matrix[x][y].v_name[i] == vertex_matrix[new_X][new_Y].v_name[i])      continue;
                else if(vertex_matrix[x][y].v_name[i] == '_')                                   continue;
                else if(vertex_matrix[new_X][new_Y].v_name[i] == '_')                           continue;
                else{   seq_eq = false;     break;}
            }


            if(seq_eq){
                if(used[x] == 0){
                    vertex_matrix[x][y].v_uni_code = biggest;
                    vertex_matrix[new_X][new_Y].v_uni_code = biggest;

                    cout<<"X = "<<x<<", Y = "<<y<<"\t\t";
                    for(int i = 0; i < window; i++) {
                        cout << vertex_matrix[x][y].v_name[i] << " ";
                    }                                   cout<<endl;
                    used[x] = 1;
                }


            }


        }
    }


    cout<<endl<<endl;


}


void clique_family(){
    int tmp = 0;
    int len = 0;
    int x_pos, y_pos;
    bool family = true;


    for(int x = 0; x < seq.size(); x++) {
        for (int y = 0; y < (seq[x].size() - window + 1); y++) {
            tmp = 0;
            family = true;
            while(family) {
                if (vertex_matrix[x][y].v_uni_code + 1  == vertex_matrix[x][y + 1].v_uni_code) {
                    tmp++;
                    y++;
                }else family = false;
            }
            if(tmp > len){
                len = tmp;
                x_pos = x;
                y_pos = y-tmp;
            }
        }
    }



    cout<<"\nX = "<<x_pos<<", Y = "<<y_pos<<", len = "<<len<<endl;
    for(int i = 0; i < len; i++){
        cout<<vertex_matrix[x_pos + i][y_pos + i].v_name[0];
        if(i+1 == len){
            for(int j = 1; j<window; j++){
                cout<<vertex_matrix[x_pos + i][y_pos + i].v_name[j];
            }
        }
    }

    cout<<endl;
}




int main() {
    window = 5;    quality_lvl = 10;    deletion_nr = 1;

    cout<<"Window = "<<window<<", quality_lvl = "<<quality_lvl<<", deletion_nr = "<<deletion_nr<<endl;

    //TAB FOR FIND MOST COMMON UNI_CODE
    for(auto & i : tab_uni_code)
        for(int & j : i)
            j=0;

    data_set_nr = 4;



    vector <Vertex> tmp_vertex;

//    window_quality()
    seq_reader();
    quality_reader();

    vertex_generator();

    edge_generator();


//    print_seq_uni_code();



    return 0;
}