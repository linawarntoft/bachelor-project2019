//Bachelor project by Lina Warntoft 2018

#ifndef NODEARR_H
#define NODEARR_H

#include <iostream>
using namespace std;

static const int VARIABLES = 11; //x, y, z, vx, vy, vz, ax, ay, az, escaper_bool, index


//---------------------- NODE CLASS FOR LINKED LIST ----------------------

template <typename T>
class Node{
public:
	Node();
	Node(const T data[], Node* link);
	
	void set_data(const T data[]);
	void set_data(const T data, const int index);
	void set_link(Node* link);
	
	T get_data(const int index) const;
	Node* get_link() const;

private:
	T data_[VARIABLES];
	Node* link_;
};

//------------------------------------------------------------------------
//----------------------- ITERATOR FOR LINKED LIST -----------------------

// It turned out not to be needed
/*
class Iterator{
public:
	Iterator(Node initial=nullptr);
	double operator*() const;
	Iterator& operator++();
	Iterator& operator++(int);
	bool operator==(const Iterator other) const;
	bool operator!=(const Iterator other) const;
	
private:
	Node current_;
};*/
//------------------------------------------------------------------------
//------------------- FUNCTIONS RELEVANT TO LINKED LIST ------------------
template <typename T>
void list_clear(Node<T>*& head_ptr);
template <typename T>
void list_head_remove(Node<T>*& head_ptr);

//------------------------------------------------------------------------
//---------------------------- NODE FUNCTIONS ----------------------------

template <typename T>
Node<T>::Node(){
	for(int i=0; i < VARIABLES; i++){
		data_[i] = 0;
	}
	link_ = nullptr;
}

template <typename T>
Node<T>::Node(const T data[], Node<T>* link){
	for (int i=0; i < VARIABLES; i++){
		data_[i] = data[i];
	}
	link_ = link;
}

template <typename T>
void Node<T>::set_data(const T data[]){
	for (int i=0; i < VARIABLES; i++){
		data_[i] = data[i];
	}
}

template <typename T>
void Node<T>::set_data(const T data, const int index){
	data_[index] = data;
}

template <typename T>
void Node<T>::set_link(Node<T>* link){
	link_ = link;
}

template <typename T>
T Node<T>::get_data(const int index) const{
	return data_[index];
}

template <typename T>
Node<T>* Node<T>::get_link() const{
	return link_;
}

template <typename T>
void list_clear(Node<T>*& head_ptr){
    while (head_ptr != nullptr){
        list_head_remove(head_ptr);
	}
}

template <typename T>
void list_head_remove(Node<T>*& head_ptr){
    Node<T>* remove_ptr;
    
    remove_ptr = head_ptr;
    head_ptr = head_ptr->get_link();
    delete remove_ptr; //delete[] remove_ptr;
}

//------------------------------------------------------------------------
//-------------------------- ITERATOR FUNCTIONS --------------------------
/*
Iterator::Iterator(Node* initial=nullptr){
	current_ = initial;
}
	
double Iterator::operator*() const{
	return current_->get_data();
}

Iterator& Iterator::operator++(){
	current_ = current_->get_link();
	return *this;
}

Iterator& Iterator::operator++(int){
	Iterator temp = *this;
	current_ = current_->get_link();
	return temp;
}

bool Iterator::operator==(const Iterator other) const{
	return current_ == other.current_;
}

bool Iterator::operator!=(const Iterator other) const{
	return current_ != other.current_;
}*/
	
//------------------------------------------------------------------------


#endif //NODEARR_H